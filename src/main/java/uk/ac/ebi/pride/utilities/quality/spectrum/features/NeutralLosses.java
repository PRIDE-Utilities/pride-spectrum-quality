package uk.ac.ebi.pride.utilities.quality.spectrum.features;


import uk.ac.ebi.pride.utilities.data.core.Spectrum;
import uk.ac.ebi.pride.utilities.mol.NeutralLoss;
import uk.ac.ebi.pride.utilities.quality.utils.SpectrumFeatureType;

import java.util.HashMap;
import java.util.Map;

/**
 * This class looks for mass differences
 * @author Franz Roos & Jonas Grossmann
 */

public class NeutralLosses implements FeatureCalculator{


    private static NeutralLosses instance = null;

    Map<SpectrumFeatureType, Object> features = new HashMap<SpectrumFeatureType, Object>();

    protected NeutralLosses(){
        features.put(SpectrumFeatureType.QUALSCORE_AMONIA_17_SC,    0);
        features.put(SpectrumFeatureType.QUALSCORE_AMONIA_17_AVGB,  0);
        features.put(SpectrumFeatureType.QUALSCORE_AMONIA_17_AVGMB, 0);
        features.put(SpectrumFeatureType.QUALSCORE_AMONIA_18_SC,    0);
        features.put(SpectrumFeatureType.QUALSCORE_AMONIA_18_AVGB,  0);
        features.put(SpectrumFeatureType.QUALSCORE_AMONIA_18_AVGMB, 0);
        features.put(SpectrumFeatureType.QUALSCORE_AMONIA_28_SC,    0);
        features.put(SpectrumFeatureType.QUALSCORE_AMONIA_28_AVGB,  0);
        features.put(SpectrumFeatureType.QUALSCORE_AMONIA_28_AVGMB, 0);
    }

    public static NeutralLosses getInstance(){
        if(instance == null)
            instance = new NeutralLosses();
        return instance;
    }

    /**
     *
     * @param spectrum Spectrum
     * @param charge charge
     * @return
     */
    @Override
    public Map<SpectrumFeatureType, Object> computeFeature(Spectrum spectrum, int charge) {

        if (spectrum == null || spectrum.getMassIntensityMap().length == 0)
            return features;

        double[][] peakList = spectrum.getMassIntensityMap();

        int repeats = 8;
        double sum = 0, avg = 0;


        // average background neutral losses over several nonsense offsets such as -5, -10, -15, ...
        for(int i = 1; i <= repeats; i++){
            sum += lookForMassDifference(peakList, 0.5, repeats * (-5), false);
        }
        avg = sum/repeats;

        // root transformation for a more gaussian distribution
        features.put(SpectrumFeatureType.QUALSCORE_AMONIA_17_SC,    Math.sqrt(lookForMassDifference(peakList, 0.5, 17, false)));
        features.put(SpectrumFeatureType.QUALSCORE_AMONIA_17_AVGB,  Math.sqrt(lookForMassDifference(peakList, 0.5, 17, false) / (avg + 1)));
        features.put(SpectrumFeatureType.QUALSCORE_AMONIA_17_AVGMB, Math.log(Math.abs(lookForMassDifference(peakList, 0.5, 17, false) - (avg)) + 1));
        features.put(SpectrumFeatureType.QUALSCORE_AMONIA_18_SC,    Math.sqrt(lookForMassDifference(peakList, 0.5, 18, false)));
        features.put(SpectrumFeatureType.QUALSCORE_AMONIA_18_AVGB,  Math.sqrt(lookForMassDifference(peakList, 0.5, 18, false) / (avg + 1)));
        features.put(SpectrumFeatureType.QUALSCORE_AMONIA_18_AVGMB, Math.log(Math.abs(lookForMassDifference(peakList, 0.5, 18, false) - (avg)) + 1));
        features.put(SpectrumFeatureType.QUALSCORE_AMONIA_28_SC,    Math.sqrt(lookForMassDifference(peakList, 0.5, 28, false)));
        features.put(SpectrumFeatureType.QUALSCORE_AMONIA_28_AVGB,  Math.sqrt(lookForMassDifference(peakList, 0.5, 28, false) / (avg + 1)));
        features.put(SpectrumFeatureType.QUALSCORE_AMONIA_28_AVGMB, Math.log(Math.abs(lookForMassDifference(peakList, 0.5, 28, false) - (avg)) + 1));
        return features;
    }

    private static double lookForMassDifference(double[][] peakList, double masstolerance, double target_massdiff, boolean weighted) {
        double count = 0;
        int lower_compl;
        double deviation = 0;
        
        // compare all against all peaks, but leave out those that are out of the mass tolerance anyway
        lower_compl = 0;
        
        for (int i = 0; i < peakList.length; i++) {
            for (int j = lower_compl; j < peakList.length; j++) {
                deviation = peakList[j][0] - peakList[i][0] - target_massdiff;
                if(Math.abs(deviation) < masstolerance){
                    if(weighted == true){
                        count += Math.sqrt(peakList[i][1] * peakList[j][1]);
                    }else{
                        count += 1;
                    }
                }else{
                    
                    // try to save time
                    if(deviation < 0){
                        lower_compl += 1; // no point in searching as low down as here
                    }else{
                        break; // no point in searching as far up as here
                    }
                }
            }
        }
        
        return count;
    }
    
}
