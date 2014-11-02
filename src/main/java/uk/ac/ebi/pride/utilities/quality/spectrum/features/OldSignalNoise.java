
 /**
 * This class looks for real isotope patterns in a spectrum and compares them
 * to a theoretical isotope pattern calculated via a Poisson distribution.
 * Only the main peak and peaks at +1, +2, +3 Da are considered
 *
 * Franz Roos
 * August 13, 2004
 * */

package uk.ac.ebi.pride.utilities.quality.spectrum.features;


 import uk.ac.ebi.pride.utilities.data.core.Spectrum;
 import uk.ac.ebi.pride.utilities.quality.spectrum.common.SpectrumUtils;
 import uk.ac.ebi.pride.utilities.quality.utils.SpectrumFeatureType;

 import java.util.HashMap;
 import java.util.Map;

 public class OldSignalNoise implements FeatureCalculator{

     @Override
     public Map<SpectrumFeatureType, Object> computeFeature(Spectrum spectrum, int charge) {
         Map<SpectrumFeatureType, Object> features = new HashMap<SpectrumFeatureType, Object>();

         double[][] peakList = spectrum.getMassIntensityMap();
         double score1 = 0;
         double score2 = 0;
         double score3 = 0;
         double cutoff = SpectrumUtils.getTotIonCurrentCount(spectrum) / (spectrum.getMassIntensityMap().length + 1);

         for (int i = 0; i < peakList.length; i++) {
             if(peakList[i][1] > cutoff){ // signal peak
                 score1 += 1;
             }else{ // noise peak
                 score2 += 1;
             }
         }

         score3 = Math.sqrt(Math.sqrt(score1 / (score2 + 1))) * 10;  // double root transformation for distribution plot
         score1 = Math.sqrt(Math.sqrt(score1)) * 10;                 // double root transformation for distribution plot
         score2 = Math.sqrt(Math.sqrt(score2)) * 10;                 // double root transformation for distribution plot

         features.put(SpectrumFeatureType.QUALSCORE_ISOTOPE_SIGN, score1);
         features.put(SpectrumFeatureType.QUALSCORE_ISOTOPE_NOISE, score2);
         features.put(SpectrumFeatureType.QUALSCORE_ISOTOPE_STN, score3);

         return features;
     }
 }
