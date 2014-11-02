package uk.ac.ebi.pride.utilities.quality.spectrum.features;


import uk.ac.ebi.pride.utilities.data.controller.DataAccessUtilities;
import uk.ac.ebi.pride.utilities.data.core.Spectrum;
import uk.ac.ebi.pride.utilities.quality.utils.SpectrumFeatureType;

import java.util.HashMap;
import java.util.Map;

/**
 * This class generates additional scores based on the complementarity
 * between b- and y-ions.
 * It is compatible to triply charged spectra.
 * @author Franz Roos & Jonas Grossmann
 */

public class TripleChargedComplementarity implements FeatureCalculator{

    private static TripleChargedComplementarity instance = null;

    private Map<SpectrumFeatureType, Object> features = new HashMap<SpectrumFeatureType, Object>();

    protected TripleChargedComplementarity(){
        //In some cases we need to return default values.
        features.put(SpectrumFeatureType.QUALSCORE_Complement_B_Y_IONS_3CHARGE_SC, 0);
        features.put(SpectrumFeatureType.QUALSCORE_Complement_B_Y_IONS_3CHARGE_AB, 0);
        features.put(SpectrumFeatureType.QUALSCORE_Complement_B_Y_IONS_3CHARGE_SA, 0);
    }

    public static TripleChargedComplementarity getInstance(){
        if(instance == null)
            instance = new TripleChargedComplementarity();

        return instance;
    }
    
    @Override
    public Map<SpectrumFeatureType, Object> computeFeature(Spectrum spectrum, int charge) {

        if (spectrum == null || spectrum.getMassIntensityMap().length == 0)
            return features;

        double[][] peakList = spectrum.getMassIntensityMap();

        double[] count = new double[peakList.length];
        double precursor_mz = DataAccessUtilities.getPrecursorMz(spectrum);

        int repeats = 8;
        double tolerance = 0.9;
        double score1 = 0, score2 = 0, score3 = 0, background_score = 0;

        // average background complementarity count over several nonsense offsets such as -5, -10, -15, ...
        for (int i = 1; i <= repeats; i++) {
            background_score += comparePeaks(peakList, tolerance, precursor_mz + (i - repeats/2 - 0.5) * 10);
        }

        score1 = Math.log(comparePeaks(peakList, tolerance, precursor_mz + 0) + 1E-4); // simple count
        score2 = comparePeaks(peakList, tolerance, precursor_mz + 0) - (background_score / repeats); // subtract background to normalize
        score3 = Math.sqrt(comparePeaks(peakList, tolerance, precursor_mz + 0) / ((background_score / repeats) + 1) ); // divide by background to normalize

        features.put(SpectrumFeatureType.QUALSCORE_Complement_B_Y_IONS_3CHARGE_SC,score1);
        features.put(SpectrumFeatureType.QUALSCORE_Complement_B_Y_IONS_3CHARGE_AB, score3);
        features.put(SpectrumFeatureType.QUALSCORE_Complement_B_Y_IONS_3CHARGE_SA, score2);

        return features;
    }


    // method to count complementary peaks for a given parent mass; running time optimized
    private static double comparePeaks(double[][] peakList, double masstolerance, double precursor_mz) {
        int upper_compl;
        int count = 0, count_doubly = 0, count_triply = 0, count_triply21 = 0;
        double massdiff_doubly = 0, massdiff_triply = 0, massdiff_triply21 = 0;
        
        // doubly charged spectra: compare singly charged daughter ions versus singly charged daughter ions
        upper_compl = peakList.length - 1;
        for (int i = 0; i < peakList.length; i++) {
            for (int j = upper_compl; j >= 0; j--) { // going all the way down from the parent mass
                massdiff_doubly = 1 * peakList[i][0] + 1 * peakList[j][0] - precursor_mz * 2;
                if (Math.abs(massdiff_doubly) <= masstolerance) {
                    count_doubly += 1; // count as correct complement
                } else {
                    // try to save time
                    if (massdiff_doubly > masstolerance && upper_compl > j) {
                        upper_compl -= 1; // no point in searching higher than this
                    } else if (massdiff_doubly < masstolerance) {
                        j = 0; // no point in searching lower than this
                    }
                }
            }
        }
        
        // triply charged spectra
        upper_compl = peakList.length - 1;
        for (int i = 0; i < peakList.length; i++) {
            for (int j = upper_compl; j >= 0; j--) { // going all the way down from the parent mass
                massdiff_triply = 1 * peakList[i][0] + 2 * peakList[j][0] - precursor_mz * 3;
                if (Math.abs(massdiff_triply) <= masstolerance) {
                    count_triply += 1; // count as correct complement
                } else {
                    // try to save time
                    if (massdiff_triply > masstolerance && upper_compl > j) {
                        upper_compl -= 1; // no point in searching higher than this
                    } else if (massdiff_triply < masstolerance) {
                        j = 0; // no point in searching lower than this
                    }
                }
            }
        }
        
        count = count_doubly > count_triply ? count_doubly : count_triply;
        //count = count_doubly + count_triply; //use maximum?
        return count;
        
    } //end comparePeaks
    
} //end class

