package uk.ac.ebi.pride.utilities.quality.spectrum.peakselectors;

import uk.ac.ebi.pride.utilities.data.core.Spectrum;

/**
 * The simplest PeakSelector: all peaks that are input are outputted
 * too.
 * 
 * @author ypriverol
 */
public class AllPeaks implements PeakSelector{

	public Spectrum transform(Spectrum spectrum){
		return spectrum;
	}
	
	public String getDescription(){
		return "All peaks";
	}
	
	public String getLabel(){
		return "All peaks";
	}
}
