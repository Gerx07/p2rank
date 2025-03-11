package cz.siret.prank.utils

import org.apache.commons.lang3.exception.ExceptionUtils

/**
 *
 */
class ErrorUtils {

    static List<String> getAllCauseMessages(Throwable throwable) {
        List<String> messages = new ArrayList<>()
        while (throwable != null) {
            messages.add(throwable.getMessage())
            throwable = throwable.getCause()
        }
        return messages
    }

    static List<String> getAllCauseMessagesWithClasses(Throwable throwable) {
        List<String> messages = new ArrayList<>()
        while (throwable != null) {
            messages.add("$throwable.class.simpleName: $throwable.message")
            throwable = throwable.getCause()
        }
        return messages
    }

    static String getRootCauseMessage(Throwable throwable) {
        return ExceptionUtils.getRootCauseMessage(throwable)
    }

}
