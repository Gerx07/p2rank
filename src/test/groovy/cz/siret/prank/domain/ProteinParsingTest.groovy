package cz.siret.prank.domain

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import org.junit.jupiter.api.Test

/**
 *
 */
@Slf4j
@CompileStatic
class ProteinParsingTest {

    static final String DIR = 'src/test/resources/data/tricky_cases/parsing'


    @Test
    void testParsePdbeUpdatedCif() {
        Protein protein = Protein.load("$DIR/4gqq.cif")

        Protein protein2 = Protein.load("$DIR/4gqq_updated.cif")   // failed in 2.5
    }

    @Test
    void testParsePdbeUpdatedCif2() {
        Protein protein = Protein.load("$DIR/1fbl.cif")

        Protein protein2 = Protein.load("$DIR/1fbl_updated.cif")  // OK in in 2.5
    }

}
