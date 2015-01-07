package uk.ac.ebi.pride.utilities.quality.spectrum.peakselectors;


import uk.ac.ebi.pride.utilities.data.core.Spectrum;

/**
 * PeakSelector is an interface that allows a class to
 * transform a spectrum into a filtered Spectrum. This is used
 * to calculate the same set of features for different
 * aspects of spectra.
 * 
 * @author @ypriverol
 */
public interface PeakSelector{
    /**
     * Transform one spectrum in another new one by a transformation process (peak selection, etc)
     * @param spectrum
     * @return
     */
	public Spectrum transform(Spectrum spectrum) throws CloneNotSupportedException;

	/**
     * Every method has they are own label to know which method was applied.
     * @return
     */
    public String getDescription();

    /**
     * Get a code for a feature, like an identifier  for the feature.
     * @return
     */
    public String getCode();


}
