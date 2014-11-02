package uk.ac.ebi.pride.utilities.quality.utils;

/**
 * Created with IntelliJ IDEA.
 * @author yperez
 */

public class PrideFeatureException extends RuntimeException {

    public PrideFeatureException(String message) {
        super(message);
    }

    public PrideFeatureException(String message, Throwable cause) {
        super(message, cause);
    }
}