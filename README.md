# IIR_Elliptic_Filter_CS_Library

C# Library to calculate the coefficients of the Elliptic filter and to filter the data.

This code calculates the coefficients of the Band-pass, Band-stop, Low-pass and High-pass Elliptic filters. It also filters the data, but no zero-phase delay is applied. The name space is: IIR_Elliptic_Filter. The code follows the same steps as in Matlab.

Each filter function will return a 2 rows x N coefficients 2D vector, where Row 1 = Numerator and Row 2 = Denumerator. The method "check_stability_iir" can be used to check the stability of the filter.

Band-pass: the function is "double[][] Lp2bp (int, double, double, double, double)". The first three arguments are the order of the filter, the decibels of peak - to - peak passband ripple and the decibels of stopband attenuation down from the peak passband value, respectively. The last two arguments are the two normalized cut-off frequencies (f1/NF, f2/NF), where NF is the Nyquist frequency. This means that the cutoff frequencies must be within the interval of (0,1). Please, keep in mind that if you enter order_filt = 2, the order of the filter will be 2 * order_filt = 4;

Band-stop: the function is " double[][] Lp2bs (int, double, double, double, double)". The first three arguments are the order of the filter, the decibels of peak - to - peak passband ripple and the decibels of stopband attenuation down from the peak passband value, respectively. The last two arguments are the two normalized cut-off frequencies (f1/NF, f2/NF), where NF is the Nyquist frequency. This means that the cutoff frequencies must be within the interval of (0,1). Please, keep in mind that if you enter order_filt = 2, the order of the filter will be 2 * order_filt = 4;

High-pass: the function is " double[][] Lp2hp (int, double, double, double)". The first three arguments are the order of the filter, the decibels of peak - to - peak passband ripple and the decibels of stopband attenuation down from the peak passband value, respectively. The last argument is the normalized cut-off frequency (f/NF), where NF is the Nyquist frequency. This means that the cutoff frequencies must be within the interval of (0,1);

Low-pass: the function is " double[][] Lp2lp (int, double, double, double)". The first three arguments are the order of the filter, the decibels of peak - to - peak passband ripple and the decibels of stopband attenuation down from the peak passband value, respectively. The last argument is the normalized cut-off frequency (f/NF), where NF is the Nyquist frequency. This means that the cutoff frequencies must be within the interval of (0,1);

Check the stability of the filter: the method is " bool Check_stability_iir (double[][] coeff_filt)". The argument is the 2D array containing the filter coefficients. It returns "true" if the filter is stable, "false" if it is unstable.

Filter the data: the method is "double[] Filter_Data(double[][] coeff_filt, double[] pre_filt_signal)". The two arguments are the filter coefficients and the signal to be filtered. It returns the filtered signal.

If you have any question and/or want to report bugs, please e-mail me (Ale) at: pressalex@hotmail.com

