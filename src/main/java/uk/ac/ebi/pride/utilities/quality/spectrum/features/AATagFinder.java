package uk.ac.ebi.pride.utilities.quality.spectrum.features;

import uk.ac.ebi.pride.utilities.data.core.Spectrum;
import uk.ac.ebi.pride.utilities.quality.utils.Constants;
import uk.ac.ebi.pride.utilities.quality.utils.SpectrumFeatureType;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

/**
 * AAFinder is a class that calculates the features that have to do
 * with checking if distances between peaks could be amino-acids.
 *
 * <P>
 * At this moment 4 values are returned by the method <CODE>consecutiveScore</code>
 * (this is also reflected by the method <CODE>getNames</CODE>):
 * <OL>
 * <LI>Plain number of peak-distances that could be amino acids (normalized by dividing this number
 * by the square of the total number of peaks)
 * <LI>Features 2 and 3 are based on chains of consecutive amino-acid-like distances: for every
 * peak <I>i</I> the number D<sub>i</sub> will represent the maximum number of steps of
 * amino-acid distance can be taken towards lower m/z values.<BR>
 * Feature 2 is the average of all these values.
 * <LI>Feature 3 is the maximum over all <I>i</I> of D<sub>i</sub>
 * <LI>Feature 4 is roughly equivalent to feature 1 (number of AA distances) but it is an abundance-weighted version (FR)
 * </OL>
 *
 * originally developed by M. Vogelzang
 * @author ypriverol
 *
 */
public class AATagFinder implements FeatureCalculator{
    /**
     * The tolerance that is used when comparing a peak distance and a possible amino acid.
     */
    private static final double TOLERANCE = 0.7;

    private static AATagFinder instance;

    Map<SpectrumFeatureType, Object> features = new HashMap<SpectrumFeatureType, Object>();

    protected AATagFinder(){
        features.put(SpectrumFeatureType.QUALSCORE_AA_MASS_DIFF,         0);
        features.put(SpectrumFeatureType.QUALSCORE_AA_MASS_ABUNDANCE_WEIGHTED,0);
        features.put(SpectrumFeatureType.QUALSCORE_AA_MASS_TAG_AVG_LONGER, 0);
        features.put(SpectrumFeatureType.QUALSCORE_AA_MASS_TAG_LONGER, 0);

    }

    public static AATagFinder getInstance(){
        if(instance == null)
            instance = new AATagFinder();
        return instance;
    }

    /**
     * Calculate the 4 above-mentioned values for a certain Scan
     * @param spectrum The Scan to calculate the values for.
     * @param charge If this is a singly or multiply charged spectrum
     * @return The calculated features
     */

    @Override
    public Map<SpectrumFeatureType, Object> computeFeature(Spectrum spectrum, int charge) {

        if (spectrum == null || spectrum.getMassIntensityMap().length == 0)
            return features;

        double[][] peakList = spectrum.getMassIntensityMap();
        double[]   count = new double[peakList.length];
            
        ArrayList[] destinations = new ArrayList[peakList.length];
        double result1 = 0, result1_pos = 0, result1_tot = 0;
        double result2 = 0;
        double result3 = 0;
        double result4 = 0, result4_pos = 0, result4_tot = 0;

        for (int i = 0; i < peakList.length; i++) {
            destinations[i] = new ArrayList();
            int aaIndex = 0, aaIndex2 = 0;
            for (int j = i + 1; j < peakList.length; j++) {
                while (aaIndex != Constants.AAMasses.length - 1 && Math.abs(peakList[i][0] + Constants.AAMasses[aaIndex + 1] - peakList[j][0]) < Math.abs(peakList[i][0] + Constants.AAMasses[aaIndex] - peakList[j][0])) {
                    aaIndex++;
                }
                if (Math.abs(peakList[i][0] + Constants.AAMasses[aaIndex] - peakList[j][0]) < TOLERANCE) {
                    destinations[i].add(new Integer(j));
                    result1_pos += 1;
                    result4_pos += (peakList[i][1] * peakList[j][1]);
                }
                result1_tot += 1;  // count all comparisons as a reference value
                result4_tot += (peakList[i][1] * peakList[j][1]);  // sum up all comparisons as a reference value
            }
        }

        if (peakList.length != 0) {
                result1 /= peakList.length*(double)peakList.length; // normalize result1
                for (int i = 0; i < peakList.length; i++) count[i] = 1;
                result2 = 0;
                for (int j = peakList.length - 2; j >= 0; j--) {
                    for (int i = 0; i < destinations[j].size(); i++) {
                        int dest = ((Integer) destinations[j].get(i)).intValue();
                        if (count[j] < count[dest] + 1) {
                            count[j] = count[dest] + 1;
                            if (count[j] > result3) result3 = count[j];
                        }
                    }
                }
                for (int i = 0; i<peakList.length; i++) result2 += count[i];
                result2 /= peakList.length;
            }
            
            //calculate final results
            result1 = result1_pos / (result1_tot + 1) * 10;
            result4 = ( result4_pos/(result1_pos + 1) ) / ( (result4_tot / (result1_tot + 1)) + 1 ) * 10; // average of positive comparisons divided by average of total comparisons
            features.put(SpectrumFeatureType.QUALSCORE_AA_MASS_DIFF,result1);
            features.put(SpectrumFeatureType.QUALSCORE_AA_MASS_ABUNDANCE_WEIGHTED,result4);
            features.put(SpectrumFeatureType.QUALSCORE_AA_MASS_TAG_AVG_LONGER, result2);
            features.put(SpectrumFeatureType.QUALSCORE_AA_MASS_TAG_LONGER, result3);
            return features;
        }

}
