package uk.ac.ebi.pride.utilities.quality.spectrum.qualityassessment;

import uk.ac.ebi.pride.utilities.data.controller.DataAccessUtilities;
import uk.ac.ebi.pride.utilities.data.core.Spectrum;
import uk.ac.ebi.pride.utilities.quality.spectrum.common.SpectrumUtils;
import uk.ac.ebi.pride.utilities.quality.utils.math.RobustMath;

import java.util.*;

/**
 * This class contains a set of spectrum functions that can be use to handler an process
 * mass spectrum.
 *
 * User: @ypriverol
 */

public class QualityAssessment {

    protected final static int TOLERANCE = 2;

    public  final static String COMPLEMENT_SINGLE_SCORE        = "Complements: simple count";
    public  final static String COMPLEMENT_DIVIDED_AVG_SCORE   = "Complements: divide by averaged background";
    public  final static String COMPLEMENT_SUBSTRACT_AVG_SCORE = "Complements: subtract averaged background";

    /**
     * Function to compute the XArea of an Spectrum following the Nu et al approach:
     * Na S, Paek E: Quality Assessment of Tandem Mass Spectra Based on Cumulative Intensity Normalization.
     * Journal of Proteome Research 2006, 5(12):3241-3248.
     * @param spectrum
     * @return the XArea Score
     */
    public static Double computeXXArea(Spectrum spectrum){

        double[] intensityClone = spectrum.getIntensityBinaryDataArray().getDoubleArray();

        Arrays.sort(intensityClone);

        List<Double> cumulativeIntensities = new ArrayList<Double>(intensityClone.length);

        double totalIntensity = 0.0;

        for(int i =0 ; i < intensityClone.length; i++) totalIntensity += intensityClone[i];

        double cumulativeIntensity = 0;
        for(int i = 0; i < intensityClone.length; i++){
            cumulativeIntensity += intensityClone[i];
            cumulativeIntensities.add(cumulativeIntensity);
        }

        //And now compute the XArea
        double triangleArea = ((double)cumulativeIntensities.size() * (Collections.max(cumulativeIntensities) - Collections.min(cumulativeIntensities)) + Collections.max(cumulativeIntensities)) / 2;

        double cumulativeArea = 0.0;

        for(Double cumulativeAreaValue: cumulativeIntensities) cumulativeArea += cumulativeAreaValue;


        //Now lets do some normalization so that the area of the triangle is always 1;
        //Modification introduced by PCC
        cumulativeArea /= triangleArea;

        //double XArea = triangleArea - cumulativeArea;
        double XArea = cumulativeArea;

        return XArea;

    }

    /**
     * Compute the Total Intensity for an Spectrum
     * @param spectrum
     * @return
     */
    public static Double computeTIC(Spectrum spectrum){
        Double tic = 0.0;
        for(int i = 0; i < spectrum.getIntensityBinaryDataArray().getDoubleArray().length; i++)
            tic = spectrum.getIntensityBinaryDataArray().getDoubleArray()[i] + tic;
        return tic;
    }

    public static void computeComplementaryIonScore(Spectrum spectrum, double intensityPercentageCuttoff){

        double[] intensityClone = spectrum.getIntensityBinaryDataArray().getDoubleArray();

        double maxIntensity = intensityClone[0];

        for(int i = 1; i < intensityClone.length; i++){
            maxIntensity = (maxIntensity > intensityClone[i]? maxIntensity:intensityClone[i]);
        }

        List<Double> mzs = new ArrayList<Double>();
        for(int i = 0; i < intensityClone.length; i ++){
            if(intensityClone[i] > (maxIntensity * intensityPercentageCuttoff)){
                mzs.add(spectrum.getMzBinaryDataArray().getDoubleArray()[i]);
            }
        }
        /*
        ions.Sort((a, b) => a.Intensity.CompareTo(b.Intensity));

        List<int> rankScore = new List<int>();

        string [] zLine = Regex.Split(ms.ZLines[0], "\t");

        double summedRanks = 0;
        ]

        double precursorMH =
        precursorMH += 1.0072;

        for (int i = 0; i < ions.Count; i++)
        {
            summedRanks += i;
            for (int j = i + 1; j < ions.Count; j++)
            {
                if (PatternTools.pTools.PPM(ions[i].MZ + ions[j].MZ, precursorMH) < ppmTolerance)
                {
                    rankScore.Add(i);
                    rankScore.Add(j);
                }
            }
        }

        rankScore = rankScore.Distinct().ToList();

        myScores.ComplementaryIonScoe = (double)rankScore.Sum() / summedRanks;
    }                                                                         */

    }

    // main method
    public static Map<String, Double> complementScores(Spectrum spectrum, boolean singlyCharged) {

        double [][] peakList = spectrum.getMassIntensityMap();

        double[] count = new double[peakList.length];

        Integer charge = DataAccessUtilities.getPrecursorCharge(spectrum.getPrecursors());
        charge = (charge == null)? DataAccessUtilities.getPrecursorChargeParamGroup(spectrum):charge;
        Double precursorMZ = DataAccessUtilities.getPrecursorMz(spectrum);
        precursorMZ = (precursorMZ == -1)? DataAccessUtilities.getPrecursorMz(spectrum):precursorMZ;

        int repeats = 8;
        double tolerance = 0.9;
        double score1 = 0, score2 = 0, score3 = 0, background_score = 0;

        // average background complementarity count over several nonsense offsets such as -5, -10, -15, ...
        for (int i = 1; i <= repeats; i++) {
            background_score += comparePeaks(peakList, tolerance, precursorMZ + (i - repeats/2 - 0.5) * 10);
        }

        score1 = Math.log(comparePeaks(peakList, tolerance, precursorMZ + 0) + 1E-4); // simple count
        score2 = comparePeaks(peakList, tolerance, precursorMZ + 0) - (background_score / repeats); // subtract background to normalize
        score3 = Math.sqrt(comparePeaks(peakList, tolerance, precursorMZ + 0) / ((background_score / repeats) + 1) ); // divide by background to normalize

        Map<String , Double> scores = new HashMap<String, Double>();
        scores.put(COMPLEMENT_SINGLE_SCORE,score1);
        scores.put(COMPLEMENT_DIVIDED_AVG_SCORE,score3);
        scores.put(COMPLEMENT_SUBSTRACT_AVG_SCORE,score2);

        return scores;
    }

    // method to count complementary peaks for a given parent mass; running time optimized
    private static double comparePeaks(double [][] peakList, double masstolerance, double precursor_mz) {
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

    public static Double calculateComplementaryIonScore(Spectrum spectrum){

        spectrum = SpectrumUtils.sortSpectrumByMass(spectrum,false);

        double[][] peakList = spectrum.getMassIntensityMap();

        int minMz = (int)peakList[0][0];

        int maxMz = (int)(peakList[peakList.length-1][0]+1);

        boolean [] peakThere = new boolean[maxMz-minMz+1];
        for(int i = 0;i<peakList.length;i++)
            peakThere[(int)(peakList[i][0] - minMz)] = true;

        int [] dev = new int[] { -10,-15,-20,-25,-30,10,15,20,25,30};

        Integer charge = DataAccessUtilities.getPrecursorCharge(spectrum.getPrecursors());
        charge = (charge == null)? DataAccessUtilities.getPrecursorChargeParamGroup(spectrum):charge;
        Double precursorMZ = DataAccessUtilities.getPrecursorMz(spectrum);
        precursorMZ = (precursorMZ == -1)? DataAccessUtilities.getPrecursorMz(spectrum):precursorMZ;

        if(charge != null && charge == 1){
            int realHits = getSingleHits(precursorMZ,peakThere,minMz);
            Integer[] values = new Integer[dev.length];
            for(int i = 0;i<dev.length;i++)
                values[i] = new Integer(getSingleHits(precursorMZ+dev[i], peakThere, minMz));

            RobustMath.NormalDistributionParameters ndp = RobustMath.normalEstimate(values);

            double z;
            if(ndp.stddev < 0.001){
                if(realHits < ndp.mean - 0.5){
                    z = -2;
                }else{
                    if(realHits > ndp.mean +0.5){
                        z = 2;
                    }else{
                        z = 0;
                    }
                }
            }else{
                z = (realHits - ndp.mean) / ndp.stddev;
            }
            return z;

        }else if(charge != null && charge > 1){

            // Check doubly / triply charged
            int realHits = getDoubleHits(precursorMZ,peakThere,minMz);

            Integer[] values = new Integer[dev.length];
            for(int i = 0;i<dev.length;i++)
                values[i] = new Integer(getDoubleHits(precursorMZ+dev[i], peakThere, minMz));

            RobustMath.NormalDistributionParameters ndp = RobustMath.normalEstimate(values);

            double z2;
            if(realHits == 0 && ndp.mean < 0.01){
                z2 = 0;
            }else{
                if(ndp.stddev > 0.3){
                    z2 = (realHits-ndp.mean) / ndp.stddev;
                }else{
                    z2 = (realHits-ndp.mean)*3;
                }
                realHits = getTripleHits(precursorMZ, peakThere, minMz);
                for(int i = 0;i<dev.length;i++)
                    values[i] = new Integer(getTripleHits(precursorMZ+dev[i], peakThere, minMz));

                ndp = RobustMath.normalEstimate(values);
                double z3 ;
                if(realHits == 0 && ndp.mean < 0.01){
                    z3 = 0;
                }else{
                    if(ndp.stddev > 0.3){
                        z3 = (realHits-ndp.mean) / ndp.stddev;
                    }else{
                        z3 = (realHits-ndp.mean)*3;
                    }
                }

                if(z3 > z2)
                    return z3;
                else
                    return z2;
            }
        }
        return null;

    }

    protected static int getSingleHits(double target, boolean[] peakThere, int start) {
        int hits = 0;
        for(int i = 0;i<peakThere.length;i++)
            if(peakThere[i])
            {
                boolean matched = false;
                int mz = i+start;
                int j = (int)(target - TOLERANCE - mz - start + 0.5);
                if(j < 0) j =0;
                int endj = (int)(target + TOLERANCE - mz - start + 0.5);
                if(endj >= peakThere.length) endj = peakThere.length-1;
                for(;j<=endj;j++)
                    if(peakThere[j])
                        hits++;
            }
        return hits;
    }

    protected static int getDoubleHits(double target, boolean[] peakThere, int start) {
        int hits = 0;
        for(int i = 0;i<peakThere.length;i++)
            if(peakThere[i]) {
                boolean matched = false;
                double mass = (i+start)-1;
                double secondMz = ((target * 2 - 2) - mass) + 1;
                int j = (int)(secondMz - TOLERANCE - start + 0.5); // what are the 0.5?
                if(j < 0) j =0;
                int endj = (int)(secondMz + TOLERANCE - start + 0.5);
                if(endj >= peakThere.length) endj = peakThere.length-1;
                for(;j<=endj;j++)
                    if(peakThere[j])
                        matched = true;
                if(matched)
                {
                    hits ++;
                    i++;
                }
            }
        return hits;
    }

    protected static int getTripleHits(double target, boolean[] peakThere, int start) {
        int hits = 0;
        for(int i = 0;i<peakThere.length;i++)
            if(peakThere[i])
            {
                boolean matched = false;
                int mz = i+start;
                double secondMass = (target * 3 - 3) - (mz-1);
                int j = (int)((secondMass+2)/2.0 - TOLERANCE - start + 0.5);
                if(j < 0) j =0;
                int endj = (int)((secondMass+2)/2.0 + TOLERANCE - start + 0.5);
                if(endj >= peakThere.length) endj = peakThere.length-1;
                for(;j<=endj;j++)
                    if(peakThere[j])
                        matched = true;
                if(matched){
                    i++;
                    hits++;
                }

            }
        return hits;
    }

}