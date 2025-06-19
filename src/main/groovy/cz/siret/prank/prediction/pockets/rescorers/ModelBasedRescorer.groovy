package cz.siret.prank.prediction.pockets.rescorers

import cz.siret.prank.domain.Pocket
import cz.siret.prank.domain.Prediction
import cz.siret.prank.domain.labeling.LabeledPoint
import cz.siret.prank.domain.labeling.ResidueLabelings
import cz.siret.prank.features.FeatureExtractor
import cz.siret.prank.features.FeatureVector
import cz.siret.prank.features.PrankFeatureExtractor
import cz.siret.prank.features.api.ProcessedItemContext
import cz.siret.prank.geom.Atoms
import cz.siret.prank.prediction.metrics.ClassifierStats
import cz.siret.prank.prediction.pockets.PocketPredictor
import cz.siret.prank.prediction.pockets.PointScoreCalculator
import cz.siret.prank.prediction.transformation.ScoreTransformer
import cz.siret.prank.program.ml.Model
import cz.siret.prank.program.params.Parametrized
import groovyx.gpars.GParsPool
import groovy.transform.CompileStatic
import groovy.transform.TypeCheckingMode
import groovy.util.logging.Slf4j
import org.biojava.nbio.structure.Atom

import static cz.siret.prank.prediction.pockets.PointScoreCalculator.applyPointScoreThreshold

/**
 * rescorer and predictor
 * 
 * Not thread safe!
 *
 * This is the main rescore used by P2RANK to make predictions based on machine learning
 *
 */
@Slf4j
@CompileStatic
class ModelBasedRescorer extends PocketRescorer implements Parametrized  {

    private final double POSITIVE_POINT_LIGAND_DISTANCE = params.positive_point_ligand_distance

    private final PointScoreCalculator calculator = new PointScoreCalculator()

    private FeatureExtractor extractorFactory
    private Model model
    private ClassifierStats stats = new ClassifierStats()

    boolean collectPoints = params.visualizations || params.predictions
    boolean visualizeAllSurface = params.vis_all_surface

    // SAS points with ligandability score for prediction and visualization
    List<LabeledPoint> labeledPoints = new ArrayList<>()


    ModelBasedRescorer(Model model, FeatureExtractor extractorFactory) {
        this.extractorFactory = extractorFactory
        this.model = model
    }

    /**
     * @param prediction
     */
    @Override
    @CompileStatic(TypeCheckingMode.SKIP)
    void rescorePockets(Prediction prediction, ProcessedItemContext context) {

        FeatureExtractor proteinExtractor = extractorFactory.createPrototypeForProtein(prediction.protein, context)

        InstancePredictor instancePredictor = InstancePredictor.create(model, proteinExtractor)

        // PRANK (just rescoring existing pockets)
        if (!params.predictions) {
            doRescore(prediction, proteinExtractor, instancePredictor)
        }

        // compute ligandability scores of SAS points for predictions and visualization
        if (params.predictions || visualizeAllSurface) {

            FeatureExtractor extractor = (proteinExtractor as PrankFeatureExtractor).createInstanceForWholeProtein()


            int n_points = extractor.sampledPoints.points.count
            labeledPoints = new ArrayList<>(n_points)
            extractor.sampledPoints.points.each { Atom point ->
                labeledPoints.add(new LabeledPoint(point))
            }

            List<FeatureVector> vectors
            GParsPool.withPool(params.fe_threads) {
                vectors = labeledPoints.collectParallel { LabeledPoint point ->
                    extractor.calcFeatureVector(point.point)
                } as List<FeatureVector>
            }

            // classification
            double[] scores = instancePredictor.predictBatch(vectors)

            // TODO refactor: use ModelBasedPointLabeler instead of this loop
            for (int i=0; i!=n_points; ++i) {
                LabeledPoint point = labeledPoints.get(i)

                // labels and statistics
                calculator.scorePoint(point, scores[i])

                point.predicted = applyPointScoreThreshold(point.score)
                point.observed = isPositivePoint(point.point, ligandAtoms)

                if (collectStats) {
                    stats.addPrediction(point.observed, point.predicted, point.score)
                }
            }

            // generate predictions
            if (params.predictions) {
                prediction.pockets = new PocketPredictor().predictPockets(labeledPoints, prediction.protein)
                prediction.reorderedPockets = prediction.pockets
                prediction.labeledPoints = labeledPoints

                if (params.label_residues) {
                    prediction.residueLabelings = ResidueLabelings.calculate(prediction, model, extractor.sampledPoints.points, labeledPoints, context)
                }
            }
        }

        proteinExtractor.finalizeProteinPrototype()
    }

    boolean isPositivePoint(Atom point, Atoms ligandAtoms) {
        if (ligandAtoms == null || ligandAtoms.empty) {
            return false
        }
        return ligandAtoms.dist(point) <= POSITIVE_POINT_LIGAND_DISTANCE
    }

    /**
     * Rescore predictions of other methods
     * TODO refactor to use PointScoreCalculator
     */
    @CompileStatic(TypeCheckingMode.SKIP)
    private void doRescore(Prediction prediction, FeatureExtractor proteinExtractor, InstancePredictor instancePredictor) {

        proteinExtractor.prepareProteinPrototypeForPockets()

        // pocket score transformers
        ScoreTransformer probaTpTransformer = ScoreTransformer.load(params.probatp_transformer)

        for (Pocket pocket : prediction.pockets) {
            FeatureExtractor extractor = proteinExtractor.createInstanceForPocket(pocket)

            double sum = 0
            double rawSum = 0

            List<LabeledPoint> pocketLabeledPoints
            GParsPool.withPool(params.fe_threads) {
                pocketLabeledPoints = extractor.sampledPoints.points.collectParallel { Atom point ->
                    FeatureVector vector = extractor.calcFeatureVector(point)
                    double pointScore = instancePredictor.predictPositive(vector)
                    boolean predicted = applyPointScoreThreshold(pointScore)
                    boolean observed = false
                    if (collectStats) {
                        observed = isPositivePoint(point, ligandAtoms)
                    }
                    return new LabeledPoint(point, observed, predicted, pointScore)
                } as List<LabeledPoint>
            }

            for (LabeledPoint lp : pocketLabeledPoints) {
                if (collectStats) {
                    stats.addPrediction(lp.observed, lp.predicted, lp.score)
                }
                sum += calculator.transformScore(lp.score)
                rawSum += lp.score
            }

            if (collectPoints) {
                labeledPoints.addAll(pocketLabeledPoints)
            }

            double pocketScore = sum
            pocket.newScore = pocketScore
            pocket.sasPoints = extractor.sampledPoints.points
            pocket.labeledPoints = pocketLabeledPoints
            pocket.auxInfo.rawNewScore = rawSum / extractor.sampledPoints.points.count // ratio of predicted ligandable points
            pocket.auxInfo.samplePoints = extractor.sampledPoints.points.count

            if (probaTpTransformer!=null) {
                pocket.auxInfo.probaTP = probaTpTransformer.transformScore(pocketScore)
            }
        }

    }

    ClassifierStats getStats() {
        return stats
    }
    
}
