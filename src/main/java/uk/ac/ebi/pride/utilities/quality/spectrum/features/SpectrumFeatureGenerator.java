/*
 
 * Created on Feb 6, 2004
 
 */

package uk.ac.ebi.pride.utilities.quality.spectrum.features;




import uk.ac.ebi.pride.utilities.data.core.Spectrum;

import uk.ac.ebi.pride.utilities.quality.spectrum.peakselectors.AllPeaks;
import uk.ac.ebi.pride.utilities.quality.spectrum.peakselectors.IsotopePeaks;
import uk.ac.ebi.pride.utilities.quality.spectrum.peakselectors.NonNoisePeaks3;
import uk.ac.ebi.pride.utilities.quality.spectrum.peakselectors.PeakSelector;
import uk.ac.ebi.pride.utilities.quality.utils.Constants;
import uk.ac.ebi.pride.utilities.quality.utils.SpectrumFeatureType;

import java.util.HashMap;
import java.util.Map;


/**
 *
 * MasterFeatureGenerator generates all features for a certain spectrum
 *
 * such as amino acid mass differences, neutral losses, complementarity.
 *
 * From the peak list, several subsets are extracted and they are
 *
 * characterized by the ParameterGen class.
 *
 * For the amino acid mass differences, a denoised peak list subset is used.
 *
 * [FR, 29 December 2004]
 *
 *
 *
 * originally developed by M. Vogelzang
 * @author  ypriverol
 *
 */


public class SpectrumFeatureGenerator{

    private PeakSelector[] firstSelector;

    private PeakSelector[] secondSelector;

    private static SpectrumFeatureGenerator instance = null;

    private Map<Integer, Map<SpectrumFeatureType, Object>> features;

    /**
     * Protected constructor for singleton pattern without parameters. In the constructor we create the parameters to compute
     * all the features in the main feature object.
     *
     */
    protected SpectrumFeatureGenerator(){

        // [FR] selectors (peak list subsets) used for evaluation by distribution criteria
        firstSelector = new PeakSelector[]{
                new AllPeaks(),                                            // Considering all peaks
                new NonNoisePeaks3(5, 80, Constants.maxPeaksPer1000Da),    // changed to 80
                new NonNoisePeaks3(5, 70, Constants.maxPeaksPer1000Da),    // changed to 70
                new IsotopePeaks(20, false, true),                         // isotope signal, good isotope fit
                new IsotopePeaks(20, false, false),                         // isotope noise, bad isotope fit
        };

        // [FR] selectors (peak list subsets) used for sequence tags and other criteria

        secondSelector = new PeakSelector[] {
                        new NonNoisePeaks3(5, 70, Constants.maxPeaksPer1000Da), //optimized for aa Features (Selectors[0])
                        new NonNoisePeaks3(5, 70, Constants.maxPeaksPer1000Da), //optimized for ComplementFeatures (Selectors[1])
                        new NonNoisePeaks3(5, 50, Constants.maxPeaksPer1000Da)  //optimized for NeutralLosses (Selectors[2])
        };

        features = new HashMap<Integer, Map<SpectrumFeatureType, Object>>();
    }

    /**
     * Public constructor of the singleton pattern
     * @return SpectrumFeatureGenerator
     */
    public static SpectrumFeatureGenerator getInstance(){
        if(instance == null){
            instance = new SpectrumFeatureGenerator();
        }
        return instance;
    }

    /**
     * Return the number of the features computed for each Spectrum.
     * @return Feature Count for each Spectrum
     */
    public int getFeatureCount(){
        return features.size();
    }

    public Map<Integer, Map<SpectrumFeatureType, Object>> computeFeatureForSpectrum(Spectrum spectrum, int charge) throws CloneNotSupportedException {
        
        // [FR] Calculates features for standard peak list subsets
        
        Spectrum analysedSpectrum;
        
        int n = 0;
        GeneralFeatures generalFeatures = GeneralFeatures.getInstance();

        for (int i = 0; i < firstSelector.length; i++) {
            analysedSpectrum = firstSelector[i].transform(spectrum);
            Map<SpectrumFeatureType, Object> values = generalFeatures.computeFeature(analysedSpectrum, charge);
            features.put(features.size() + 1, values);
        }
        
        // define separate peak list subsets for sequence tags and other amino acid mass features

        Spectrum aaSpectrum          = secondSelector[0].transform(spectrum);		// this scan is optimized for aaFeatures
        Map<SpectrumFeatureType, Object> aaScores = AATagFinder.getInstance().computeFeature(aaSpectrum, charge);
        features.put(features.size() + 1, aaScores);


        Spectrum complementSpectrum  = secondSelector[1].transform(spectrum);		// this scan is optimized for ComplementFeatures
        Map<SpectrumFeatureType, Object> complementScores = TripleChargedComplementarity.getInstance().computeFeature(complementSpectrum, charge);
        features.put(features.size() + 1, complementScores);

        Spectrum neutralLossSpectrum = secondSelector[2].transform(spectrum);		// this scan is optimized for NeutralLoss Features
        Map<SpectrumFeatureType, Object> neutrallossescores = NeutralLosses.getInstance().computeFeature(neutralLossSpectrum, charge);
        features.put(features.size() + 1, neutrallossescores);

        return features;
    }
    
}

