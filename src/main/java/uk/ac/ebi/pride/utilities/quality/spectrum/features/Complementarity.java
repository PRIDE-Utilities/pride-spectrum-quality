package uk.ac.ebi.pride.utilities.quality.spectrum.features;


import uk.ac.ebi.pride.utilities.data.controller.DataAccessUtilities;
import uk.ac.ebi.pride.utilities.data.core.Spectrum;
import uk.ac.ebi.pride.utilities.quality.utils.SpectrumFeatureType;

import java.util.HashMap;
import java.util.Map;

/**
 * This class generates additional scores based on the complementarity
 * between b- and y-ions.
 * originally developed by Franz Roos & Jonas Grossmann
 *
 * @author ypriverol
 */

public class Complementarity implements FeatureCalculator{


    private Map<SpectrumFeatureType, Object> features = new HashMap<SpectrumFeatureType, Object>();

    private static Complementarity instance = null;

    protected Complementarity(){
        features.put(SpectrumFeatureType.QUALSCORE_Complement_B_Y_IONS_3vs3, 0);
        features.put(SpectrumFeatureType.QUALSCORE_Complement_B_Y_IONS_1vs3, 0);
        features.put(SpectrumFeatureType.QUALSCORE_Complement_B_Y_IONS_iso1, 0);
        features.put(SpectrumFeatureType.QUALSCORE_Complement_B_Y_IONS_sgn,  0);

    }

    public static Complementarity getInstance(){
        if(instance == null)
            instance = new Complementarity();
        return instance;
    }

    @Override
    public Map<SpectrumFeatureType, Object> computeFeature(Spectrum spectrum, int charge) {
        double[][] peakList = spectrum.getMassIntensityMap();
        double[] count = new double[peakList.length];
        double chargestate = 2;
        double parentmass = DataAccessUtilities.getPrecursorMz(spectrum) * chargestate;

        int repeats = 8;
        double score1 = 0, sum1 = 0;
        double score2 = 0, sum2 = 0;
        double score3 = 0, sum3 = 0;
        double score4 = 0, sum4 = 0;
        double score5 = 0;
        double score6 = 0;
        double score7 = 0;

        // average background complementarity count over several nonsense offsets such as -5, -10, -15, ...
        for(int i = 1; i <= repeats; i++){
            sum1 += comparePeaks(peakList, 1.0, parentmass - i * 5);
            sum2 += comparePeaks(peakList, 1.0, parentmass - i * 5);
            sum3 += comparePeaks(peakList, 1.0, parentmass - i * 5);
            sum4 += comparePeaks(peakList, 1.0, parentmass - i * 5);
        }

        score1 = Math.sqrt(comparePeaks(peakList, 0.3, parentmass) / ((sum1 / repeats) + 1) * 10); // avoid division by zero
        score2 = Math.sqrt(comparePeaks(peakList, 1.0, parentmass) / ((sum2 / repeats) + 1) * 10); // use -10 as background noise measurement
        score3 = Math.sqrt(comparePeaks(peakList, 0.5, parentmass + 1 ) / ((sum3 / repeats) + 1)); // isotope
        score4 = Math.sqrt(comparePeaks(peakList, 0.5, parentmass + 0 ) / ((sum4 / repeats) + 1)); // main signal
        score5 = Math.sqrt(comparePeaks(peakList, 0.5, parentmass - 17) / (comparePeaks(peakList, 0.5, parentmass - 10) + 1)); // ammonia
        score6 = Math.sqrt(comparePeaks(peakList, 0.5, parentmass - 18) / (comparePeaks(peakList, 0.5, parentmass - 10) + 1)); // water
        score7 = Math.sqrt(comparePeaks(peakList, 0.5, parentmass - 28) / (comparePeaks(peakList, 0.5, parentmass - 10) + 1)); // carbon monoxide

        Map<SpectrumFeatureType, Object> features = new HashMap<SpectrumFeatureType, Object>();
        features.put(SpectrumFeatureType.QUALSCORE_Complement_B_Y_IONS_3vs3, score1);
        features.put(SpectrumFeatureType.QUALSCORE_Complement_B_Y_IONS_1vs3, score2);
        features.put(SpectrumFeatureType.QUALSCORE_Complement_B_Y_IONS_iso1, score3);
        features.put(SpectrumFeatureType.QUALSCORE_Complement_B_Y_IONS_sgn, score4);
        return features;
//        return new double[] {score1, score2, score3, score4};

    }

     // method to count complementary peaks for a given parent mass; running time optimized
     private static double comparePeaks(double[][] peakList, double masstolerance, double parentmass){
         double count = 0;
         int upper_compl;
         double massdiff = 0;

         // compare all against all peaks, but leave out those that are out of the mass tolerance anyway
         upper_compl = peakList.length - 1;
         for (int i = 0; i < peakList.length; i++){
             for (int j = upper_compl; j > 0; j--){
                 // going all the way down from the parent mass
                 massdiff = peakList[i][0] + peakList[j][0] - parentmass;
                 if(Math.abs(massdiff) < masstolerance){
                     count += 1;  // correct complement
                 }else{
                      // try to save time
                      if(massdiff > 0){
                          upper_compl = upper_compl - 1; // no point in searching as far up as here
                      }else{
                          break; // no point in searching as low as here
                      }
                 }
             }
         }
         return count;
     }


}
