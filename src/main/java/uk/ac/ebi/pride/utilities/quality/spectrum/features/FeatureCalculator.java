package uk.ac.ebi.pride.utilities.quality.spectrum.features;

import uk.ac.ebi.pride.utilities.data.core.Spectrum;
import uk.ac.ebi.pride.utilities.quality.utils.SpectrumFeatureType;

import java.util.Map;

/**
 * @author ypriverol
 */
public interface FeatureCalculator {
    /**
     * Compute the features for an spectrum knowing the charge of it. This Interface allow to define a common structure to retrieve
     * information from spectra.
     * @param spectrum Spectrum
     * @param charge charge
     * @return Map with SpectrumFeatureType and the corresponding value it can be double or integer or a List
     */
    public Map<SpectrumFeatureType, Object> computeFeature(Spectrum spectrum, int charge);
}
