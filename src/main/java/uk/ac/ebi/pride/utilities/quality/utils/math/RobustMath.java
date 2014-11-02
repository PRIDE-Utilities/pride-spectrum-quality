/*
 * Created on Oct 29, 2003
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */
package uk.ac.ebi.pride.utilities.quality.utils.math;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

/**
 * <code>RobustMath</code> implements some methods that try to derive some
 * statistics in a robust (ie. unmodified by some sporadic outliers) way.
 *
 * @author M. Vogelzang
 */
public class RobustMath
{
	public static class NormalDistributionParameters
	{
		public double mean, stddev, variance, skewness, kurtosis;
	}

	/**
	 * Calculate median and standard deviation from the median in a sorted array of values.
	 * @param values A sorted array of values.
	 * @return A NDP with the median in the mean field and standard deviation
	 * [FR 8 Sep 2004] added skewness and kurtosis as parameters
	 */
	public static NormalDistributionParameters estimateWithMedian(Number[] values)
	{
		NormalDistributionParameters ndp = new NormalDistributionParameters();
		if(values.length == 0)
			return null;
		if(values.length % 2 == 1)
			ndp.mean = values[values.length/2].doubleValue();
		else
			ndp.mean = (values[values.length/2-1].doubleValue() + values[values.length/2].doubleValue())/2;

		double sum  =0;
		for(int i =0;i<values.length;i++)
		{
			sum += (values[i].doubleValue() - ndp.mean) * (values[i].doubleValue() - ndp.mean);
		}
		ndp.stddev = Math.sqrt(sum / values.length);
		return ndp;
	}

	public static NormalDistributionParameters estimatePercentiles(Number[] values, double percentile, int maxPeaks)
	{
		Number[] values2 = new Number[values.length];
		int index;
		for (int i = 0; i < values.length; i++){values2[i] = values[i];} // copy values to new array
		Arrays.sort(values2); // sort array to find out percentiles

		NormalDistributionParameters ndp = new NormalDistributionParameters();
		if(values.length == 0)
			return null;

		index = (int)(values2.length * percentile / 100);
		if(index < values2.length - maxPeaks) index = values2.length - maxPeaks;
		if(index < 0) index = 0;
		if(index > values2.length - 1) index = values2.length - 1;

		ndp.mean = values2[index].doubleValue(); //mean is actually not an adequate description here
		ndp.stddev = 0;
		return ndp;
	}

	public static NormalDistributionParameters normalEstimate(Number[] values)
	{
		NormalDistributionParameters ndp = new NormalDistributionParameters();
		double sum1 = 0, sum2 = 0, sum3 = 0, sum4 = 0;
		double n = values.length;
		for (int i = 0; i < n; i++)
			sum1 += values[i].doubleValue();
		ndp.mean = sum1 / n;

		for (int i = 0; i < n; i++){
			sum2 += Math.pow(values[i].doubleValue() - ndp.mean, 2);
		}
		ndp.variance = sum2 / n;
		ndp.stddev   = Math.sqrt(ndp.variance);

		if(n > 5 && ndp.stddev > 0) {  // values below may cause crashes in LDA
                  for (int i = 0; i < n; i++){
                    sum3 += Math.pow(values[i].doubleValue() - ndp.mean, 3); // for skewness
                    sum4 += Math.pow(values[i].doubleValue() - ndp.mean, 4); // for kurtosis
                  }
                  ndp.skewness = sum3 / (n * Math.pow(ndp.stddev, 3));
                  ndp.kurtosis = sum4 / (n * Math.pow(ndp.stddev, 4));
		}else{
                  ndp.skewness = 0;
                  ndp.kurtosis = 0;
		}

		return ndp;
	}

	public static NormalDistributionParameters normalEstimate(double[] values)
	{
		NormalDistributionParameters ndp = new NormalDistributionParameters();
		double sum1 = 0, sum2 = 0, sum3 = 0, sum4 = 0;
		double n = values.length;
		for (int i = 0; i < n; i++)
			sum1 += values[i];
		ndp.mean = sum1 / n;

		for (int i = 0; i < n; i++){
			sum2 += Math.pow(values[i] - ndp.mean, 2);
		}
		ndp.variance = sum2 / n;
		ndp.stddev   = Math.sqrt(ndp.variance);

		if(n > 5 && ndp.stddev > 0) {  // values below may cause crashes in LDA
                  for (int i = 0; i < n; i++){
                    sum3 += Math.pow(values[i] - ndp.mean, 3); // for skewness
                    sum4 += Math.pow(values[i] - ndp.mean, 4); // for kurtosis
                  }
                  ndp.skewness = sum3 / (n * Math.pow(ndp.stddev, 3));
                  ndp.kurtosis = sum4 / (n * Math.pow(ndp.stddev, 4));
		}else{
                  ndp.skewness = 0;
                  ndp.kurtosis = 0;
		}

		return ndp;
	}
        

	public static NormalDistributionParameters robustEstimate(
		Number[] values,
		int maxDeleteCount,
		double zLimit)
	{
		ArrayList list = new ArrayList();
		list.addAll(Arrays.asList(values));
		Collections.sort(list);
		NormalDistributionParameters ndp;

		boolean iterate = true;
		//ndp = normalEstimate((Number[]) list.toArray(new Number[0]));
		// NEW 1-19-2004: start with median
		ndp = estimateWithMedian((Number[])list.toArray(new Number[0]));
		while (iterate)
		{
			int removeCount = 0;

			// remove peaks smaller than lower limit
			while (list.size() > 0
				&& ((Number) list.get(0)).doubleValue()
					< ndp.mean - zLimit * ndp.stddev)
			{
				list.remove(0);
				removeCount++;
			}

			// remove peaks bigger than upper limit
			while (list.size() > 0
				&& ((Number) list.get(list.size() - 1)).doubleValue()
					> ndp.mean + zLimit * ndp.stddev)
			{
				list.remove(list.size() - 1);
				removeCount++;
			}

// the comparison below should be "<" instead of ">", but this leads to an endless loop
			if (removeCount > maxDeleteCount)
				iterate = true;
			else
				iterate = false;
			ndp = normalEstimate((Number[]) list.toArray(new Number[0]));
		}
		return ndp;
	}
}
