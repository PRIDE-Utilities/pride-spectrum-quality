package uk.ac.ebi.pride.utilities.quality.spectrum.peakselectors;

import uk.ac.ebi.pride.utilities.data.core.Spectrum;

import java.util.ArrayList;

/**
 * Scores each peak in a spectrum according to its isotope pattern
 * Looks for real isotope patterns in a spectrum and compares them
 * to a theoretical isotope pattern calculated via a peak mass dependent
 * Poisson distribution.
 * The current peak and peaks at +1, +2, +3 Da are considered
 * @author ypriverol
 */

public class IsotopePeaks implements PeakSelector{
    
    protected double similarityScoreCutoff;

    protected boolean continuousScore;

    protected boolean signal;
    
    public IsotopePeaks(double similarityScoreCutoff, boolean continuousScore, boolean signal) {
        this.similarityScoreCutoff = similarityScoreCutoff;
        this.continuousScore = continuousScore;
        this.signal = signal;
    }
    
    public Spectrum transform(Spectrum spectrum) throws CloneNotSupportedException {

        Spectrum analysisSpectrum = new Spectrum(spectrum);

        double[][] peakList = analysisSpectrum.getMassIntensityMap();

        ArrayList numbers = new ArrayList();
        
        double similarityScore = 0;
        double ktemp = 0;
        int    k = 0; // index for real distribution
        double currentmass1 = 0;
        double currentmass2 = 0;
        double[] idealDistribution = new double[4];   // Theoretical Poisson distribution for peaks 0..3
        double[] realDistribution  = new double[4];  // Real Poisson distribution for peaks 0..3
        
        double[][] newPeakList = new double[peakList.length][2];
        
        for (int i = 0; i < peakList.length; i++) {
            currentmass1 = peakList[i][0];
            
            for (int xid = 0; xid < 3; xid++) {
                idealDistribution[xid] = computePoisson(xid, currentmass1);
                realDistribution[xid] = 0; // initialize realDistribution
            }
            
            for (int j = i; j < peakList.length; j++) {
                currentmass2 = peakList[j][0];
                
                // calculate in what bin the peak belongs
                ktemp = Math.round(currentmass2 - currentmass1);
                if (ktemp > 3)break;
                
                k = 0;
                while (ktemp > k) {
                    k += 1;
                }
                
                // put the peak into the bin
                if (k > 3)break;
                realDistribution[k] = realDistribution[k] + peakList[j][1];
            }
            
            // compare real distribution with ideal distribution
            similarityScore = compareDistributions(idealDistribution, realDistribution);
            //System.out.println(similarityScore + "   " + currentmass1 + "   " + peakList[i][1]);
            newPeakList[i][0] = currentmass1;     // keep all peaks, transform abundances
            newPeakList[i][1] = similarityScore;  // keep all peaks, transform abundances
            if (signal == false) {
                if (similarityScore > similarityScoreCutoff) numbers.add(new Integer(i)); // only keep a subset of peaks
            }else{
                if (similarityScore < similarityScoreCutoff) numbers.add(new Integer(i)); // only keep a subset of peaks
            }
        }
        
        // add subset to new peak list
        double[][] newPeakList2 = new double[numbers.size()][2];
        for (int i = 0; i < numbers.size(); i++) {
            newPeakList2[i] = peakList[ ( (Integer) numbers.get(i)).intValue()];
        }
        
        if(continuousScore == true){
            analysisSpectrum.setMassIntensityMap(newPeakList);  // return all peaks, transformed abundances
        }else{
            analysisSpectrum.setMassIntensityMap(newPeakList2); // return subset of peaks
        }
        
        return spectrum;
    }
    
    @Override
    public String getDescription(){
        String returnstring = "";
        if(signal == true) returnstring = "Isotope signal, peaks which fit a mass-dependent isotope pattern predicted by a Poisson model";
        if(signal == false) returnstring = "Isotope noise, peaks which fit a mass-dependent isotope pattern predicted by a Poisson model";
        return returnstring;
    }

    @Override
    public String getLabel() {
        String returnstring = "";
        if(signal == true) returnstring = "Isotope signal (cutoff " + similarityScoreCutoff + ")";
        if(signal == false) returnstring = "Isotope noise, (cutoff: " + similarityScoreCutoff + ")";
        return returnstring;
    }
    
  /*
   *  Poisson model for isotopes
   *  This method calculates the theoretical relative height of an isotope peak within an isotope pattern
   *  normalization: the abundance of all the isotope peaks in the pattern adds up to 1
   *  formula: p(x, lambda) = e^(-lambda)*lambda^x/x!   (exclamation mark means factorial here)
   */
    
    
    public static float computePoisson(int x, double mass) {
        int factorial[] = {1, 1, 2, 6};      // factorials pre-calculated for fact(0)..fact(3)
        double excess = 0.6760/1000;         // this amount of heavy isotopes occurs per dalton; per 1000 Dalton, a peptide will on average contain 1.2582 H2 or C13 or O18 or S34, ...
        double lambda = mass * excess;       // lambda as in the Poisson formula
        double density = 0;                  // result, i.e. relative frequency 0..1 for an isotope peak x (e.g. 0.7 for x=0 which is the main peak)
        float density2 = 0;                  // to return a float
        
        density = (float)Math.exp(-lambda) * Math.pow(lambda, x) / factorial[x];
        density2 = (float)density;
        return density2;
    }
    
    public static float compareDistributions(double[] idealDistribution, double[] realDistribution) {
        
        double score1 = 0;
        float realSum = 0;
        double t_chisquare = 0;                   // root(square of differences between distributions)
        
        // count the content of the bins to normalize afterwards
        for (int xid = 0; xid < realDistribution.length - 1; xid++) {
            realSum += realDistribution[xid]; // sum up real distribution
        }
        
        // normalize real distribution to 1
        for (int xid = 0; xid < realDistribution.length - 1; xid++) {
            realDistribution[xid] = realDistribution[xid] / realSum;
        }
        realSum = 0;
        
        // compare ideal and real distributions, calculate squared differences
        for (int xid = 0; xid < realDistribution.length - 1; xid++) {
            t_chisquare += (Math.pow(idealDistribution[xid] - realDistribution[xid], 2) / (realDistribution[xid] + 0.001)); // square differences
        }
        score1 = t_chisquare; // t_chisquare is the smaller the better
        return (float)score1;
    }
    
}
