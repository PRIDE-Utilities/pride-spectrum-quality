package uk.ac.ebi.pride.utilities.quality.spectrum.features;

import uk.ac.ebi.pride.utilities.data.core.Spectrum;
import uk.ac.ebi.pride.utilities.quality.utils.SpectrumFeatureType;

import java.text.DecimalFormat;
import java.util.*;

/**
 * Function to compute the XArea of an Spectrum following the Nu et al approach:
 * Na S, Paek E: Quality Assessment of Tandem Mass Spectra Based on Cumulative Intensity Normalization.
 * Journal of Proteome Research 2006, 5(12):3241-3248.
 *
 * @author ypriverol
 */
public class XXArea implements FeatureCalculator {

    private static XXArea instance = null;

    private Map<SpectrumFeatureType, Object> features = new HashMap<SpectrumFeatureType, Object>();

    protected XXArea(){
        features.put(SpectrumFeatureType.XXArea,0);
    }

    public static XXArea getInstance(){
        if(instance == null)
            instance = new XXArea();
        return instance;
  }

    @Override
    public Map<SpectrumFeatureType, Object> computeFeature(Spectrum spectrum, int charge) {

        if (spectrum == null || spectrum.getMassIntensityMap().length == 0)
            return features;

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

        DecimalFormat df = new DecimalFormat("#.##");

        XArea = Double.parseDouble(df.format(XArea));

        features.put(SpectrumFeatureType.XXArea, XArea);

        return features;
    }
}
