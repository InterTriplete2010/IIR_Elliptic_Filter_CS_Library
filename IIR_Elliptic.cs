using System;
using System.Numerics;
using System.Collections.Generic;
using System.Text;
using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Factorization;
using System.Net.Http.Headers;
using System.Linq;


namespace IIR_Elliptic_Filter
{
    public class IIR_Elliptic
    {

        //Global variables
        double fs = 2;
        double u_f1;
        double u_f2;
        double Wn;
        double Bw;

        double[][] a_matlab;
        double[] b_matlab;
        double[] c_matlab;

        double d_matlab;
        double[][] t_matlab;

        Complex complex_real = new Complex(1.0, 0.0);
        Complex complex_real_2 = new Complex(2.0, 0.0);
        Complex complex_imag = new Complex(0.0, 1.0);
        Complex complex_neg_imag = new Complex(0.0, -1.0);

        //Output of the Ellipap unction
        Complex[] z_matlab_ellipap;
        Complex[] p_matlab_ellipap;
        Complex k_matlab_ellipap;

        //Output of the Bilinear transformation
        Matrix<double> t1_arma;
        Matrix<double> t2_arma;
        Matrix<double> ad_arma;
        Matrix<double> bd_arma;
        Matrix<double> cd_arma;
        Matrix<double> dd_arma;

        double[] num_filt;   //Vector where to temporarily save the numerator
        double[] den_filt;   // Vector where to temporarily save the denumerator
        double[][] save_filt_coeff;  //Matrix where to save the numerator and denominator. First row is the numerator; second row is the denominator

        double tol_matlab = 2.220446049250313E-16;  //This value is the same used by Matlab

        //Step 1: get analog, pre - warped frequencies
        public void Freq_pre_wrapped(int type_filt, double Wnf_1, double Wnf_2)
        {

            Bw = 0;

            switch (type_filt)
            {

                //Band-pass
                case 0:

                    u_f1 = 2 * fs * Math.Tan(Math.PI * Wnf_1 / fs);
                    u_f2 = 2 * fs * Math.Tan(Math.PI * Wnf_2 / fs);

                    break;

                //Band-stop
                case 1:

                    u_f1 = 2 * fs * Math.Tan(Math.PI * Wnf_1 / fs);
                    u_f2 = 2 * fs * Math.Tan(Math.PI * Wnf_2 / fs);


                    break;

                //High-pass
                case 2:

                    u_f1 = 2 * fs * Math.Tan(Math.PI * Wnf_1 / fs);

                    break;

                //Low-pass
                case 3:

                    u_f2 = 2 * fs * Math.Tan(Math.PI * Wnf_2 / fs);

                    break;

            }

        }

        //Step 2: convert to low-pass prototype estimate
        public void Wn_f1_Wn_f2(int type_filt, double u_f1, double u_f2)
        {

            switch (type_filt)
            {

                //Band-pass
                case 0:
                    Bw = u_f2 - u_f1;
                    Wn = Math.Sqrt(u_f1 * u_f2);

                    break;

                //Band-stop
                case 1:

                    Bw = u_f2 - u_f1;
                    Wn = Math.Sqrt(u_f1 * u_f2);

                    break;

                //High-pass
                case 2:

                    Wn = u_f1;

                    break;

                //Low-pass
                case 3:

                    Wn = u_f2;

                    break;

            }
        }

        //Get Landen vector of descending moduli (method returns a vector of double)
        public List<double> Landen_vector(double k_matlab)
        {

            List<double> v_matlab = new List<double>();

            if (tol_matlab < 1)
            {

                while (k_matlab > tol_matlab)
                {

                    k_matlab = Math.Pow((k_matlab / (1 + Math.Sqrt(1 - Math.Pow(k_matlab, 2)))), 2);
                    v_matlab.Add(k_matlab);

                }

            }

            else
            {

                double M_Matlab = tol_matlab;

                for (int kk = 0; kk < M_Matlab; kk++)
                {

                    k_matlab = Math.Pow((k_matlab / (1 + Math.Sqrt(1 - Math.Pow(k_matlab, 2)))), 2);
                    v_matlab.Add(k_matlab);

                }

            }

            return v_matlab;

        }

        //Get Landen vector of descending moduli (method returns a vector of complex numbers)
        public List<Complex> Landen_vector_complex(double k_matlab)
        {

            List<Complex> v_matlab = new List<Complex>();

            if (tol_matlab < 1)
            {

                while (k_matlab > tol_matlab)
                {

                    k_matlab = Math.Pow((k_matlab / (1 + Math.Sqrt(1 - Math.Pow(k_matlab, 2)))), 2);
                    v_matlab.Add(k_matlab);

                }

            }

            else
            {

                double M_Matlab = tol_matlab;

                for (int kk = 0; kk < M_Matlab; kk++)
                {

                    k_matlab = Math.Pow((k_matlab / (1 + Math.Sqrt(1 - Math.Pow(k_matlab, 2)))), 2);
                    v_matlab.Add(k_matlab);

                }

            }

            return v_matlab;

        }


        //Method required to complete Step 3
        public double Sne(double[] u_matlab, double k_matlab)
        {

            //Get Landen vector of descending moduli
            List<double> v_matlab = Landen_vector(k_matlab);


            double[] w_matlab = new double[u_matlab.Length];

            for (int kk = 0; kk < u_matlab.Length; kk++)
            {

                w_matlab[kk] = Math.Sin(u_matlab[kk] * Math.PI / 2);

            }

            //Ascending Landen / Gauss transformation
            for (int kk = v_matlab.Count - 1; kk >= 0; kk--)
            {
                for (int ll = 0; ll < w_matlab.Length; ll++)
                {

                    w_matlab[ll] = (1 + v_matlab[kk]) * w_matlab[ll] / (1 + v_matlab[kk] * Math.Pow(w_matlab[ll], 2.0));

                }
            }

            //Clean the vector that it is longer needed
            v_matlab.Clear();

            //Calculate the product of the arry elements in w_matlab
            double prod_matlab = 1;

            for (int kk = 0; kk < w_matlab.Length; kk++)
            {

                prod_matlab *= w_matlab[kk];

            }

            return Math.Pow(prod_matlab, 4.0);

        }

        //Method required to complete Step 3 for complex numbers
        public Complex Sne(Complex u_matlab, double k_matlab)
        {

            //Get Landen vector of descending moduli
            List<Complex> v_matlab = Landen_vector_complex(k_matlab);

            Complex w_matlab = Complex.Sin(u_matlab * Math.PI / complex_real_2);

            //Ascending Landen / Gauss transformation
            for (int kk = v_matlab.Count - 1; kk >= 0; kk--)
            {

                w_matlab = (complex_real + v_matlab[kk]) * w_matlab / (complex_real + v_matlab[kk] * Complex.Pow(w_matlab, 2.0));

            }

            //Clean the vector that it is longer needed
            v_matlab.Clear();

            return w_matlab;

        }


        //Method required to complete Step 3
        public double[] Cde(double[] u_matlab, double k_matlab)
        {

            //Get Landen vector of descending moduli
            List<double> v_matlab = Landen_vector(k_matlab);

            double[] w_matlab = new double[u_matlab.Length];

            for (int kk = 0; kk < u_matlab.Length; kk++)
            {

                w_matlab[kk] = Math.Cos(u_matlab[kk] * Math.PI / 2);

            }

            //Ascending Landen / Gauss transformation
            for (int kk = v_matlab.Count - 1; kk >= 0; kk--)
            {
                for (int ll = 0; ll < w_matlab.Length; ll++)
                {

                    w_matlab[ll] = (1 + v_matlab[kk]) * w_matlab[ll] / (1 + v_matlab[kk] * Math.Pow(w_matlab[ll], 2.0));

                }
            }

            return w_matlab;

        }

        //Method required to complete Step 3 (for complex numbers)
        public Complex[] Cde(Complex[] u_matlab, double k_matlab)
        {

            //Get Landen vector of descending moduli
            List<Complex> v_matlab = Landen_vector_complex(k_matlab);

            Complex[] w_matlab = new Complex[u_matlab.Length];

            for (int kk = 0; kk < u_matlab.Length; kk++)
            {

                w_matlab[kk] = Complex.Cos(u_matlab[kk] * Math.PI / (complex_real_2));

            }

            //Ascending Landen / Gauss transformation
            for (int kk = v_matlab.Count - 1; kk >= 0; kk--)
            {
                for (int ll = 0; ll < w_matlab.Length; ll++)
                {

                    w_matlab[ll] = (complex_real + v_matlab[kk]) * w_matlab[ll] / (complex_real + v_matlab[kk] * Complex.Pow(w_matlab[ll], 2.0));

                }
            }

            return w_matlab;

        }


        //Method required to complete Step 3
        public Complex Asne(Complex w_matlab, double k_matlab)
        {

            //Get Landen vector of descending moduli
            List<Complex> v_matlab = Landen_vector_complex(k_matlab);
            Complex v1_matlab;

            for (int kk = 0; kk < v_matlab.Count; kk++)
            {

                if (kk == 0)
                {

                    v1_matlab = k_matlab;

                }

                else
                {

                    v1_matlab = v_matlab[kk - 1];

                }

                w_matlab = w_matlab / (complex_real + Complex.Sqrt(complex_real - Complex.Pow(w_matlab, 2.0) * Complex.Pow(v1_matlab, 2.0))) * complex_real_2 / (complex_real + v_matlab[kk]);

            }

            Complex u_matlab = (complex_real_2 * Complex.Acos(w_matlab)) / Math.PI;

            if (u_matlab.Real == 1 && u_matlab.Imaginary == 0)
            {

                u_matlab = 0;

            }

            //Clear some memory
            v_matlab.Clear();

            double[] K_Kprime = Ellipk(k_matlab);

            double R_matlab = K_Kprime[1] / K_Kprime[0];

            double temp_Z_I = u_matlab.Real % 4.0;

            //Extract the sign
            double sign = 1;

            if (temp_Z_I < 0)
            {

                sign = -1;

            }

            //Check the absolute value
            double abs_val_check = 0;
            if (Math.Abs(temp_Z_I) > 2)
            {

                abs_val_check = 1;

            }

            double Z_I_D = temp_Z_I - temp_Z_I * sign * abs_val_check;

            double temp_Z_II = u_matlab.Imaginary % 4.0;

            //Extract the sign
            sign = 1;

            if (temp_Z_II < 0)
            {

                sign = -1;

            }

            //Check the absolute value
            abs_val_check = 0;
            if (Math.Abs(temp_Z_II) > 2 * R_matlab)
            {

                abs_val_check = 1;

            }

            double Z_II_D = temp_Z_II - temp_Z_II * sign * abs_val_check;

            Complex u_matlab_output = new Complex(Z_I_D, Z_II_D);

            return u_matlab_output;

        }

        //Complete elliptic integral of first kind
        public double[] Ellipk(double k_matlab)
        {
            double[] K_Kprime = new double[2] { 1, 1 };

            double k_min = 1E-6;
            double k_max = Math.Sqrt(1 - Math.Pow(k_min, 2));
            double kp = 0;
            double L_matlab;
            List<double> v_matlab;
            List<double> vp_matlab;

            if (k_matlab == 1)
            {

                K_Kprime[0] = double.PositiveInfinity;

            }

            else if (k_matlab > k_max)
            {

                kp = Math.Sqrt(1 - Math.Pow(k_matlab, 2));
                L_matlab = -Math.Log(kp / 4);
                K_Kprime[0] = L_matlab + (L_matlab - 1) * Math.Pow(kp, 2) / 4;

            }

            else
            {

                v_matlab = Landen_vector(k_matlab);

                for (int kk = 0; kk < v_matlab.Count; kk++)
                {

                    K_Kprime[0] *= 1 + v_matlab[kk];

                }

                K_Kprime[0] = K_Kprime[0] * (Math.PI / 2);

            }

            if (k_matlab == 0)
            {

                K_Kprime[1] = double.PositiveInfinity;

            }

            else if (k_matlab < k_min)
            {

                L_matlab = -Math.Log(k_matlab / 4);
                K_Kprime[1] = L_matlab + (L_matlab - 1) * Math.Pow(kp, 2.0) / 4;

            }

            else
            {

                kp = Math.Sqrt(1 - Math.Pow(k_matlab, 2.0));
                vp_matlab = Landen_vector(kp);

                for (int kk = 0; kk < vp_matlab.Count; kk++)
                {

                    K_Kprime[1] *= 1 + vp_matlab[kk];

                }

                K_Kprime[1] = K_Kprime[1] * Math.PI / 2;

            }

            //Clear some memory
            //v_matlab.Clear();
            //vp_matlab.Clear();

            return K_Kprime;

        }

        //Sort complex numbers into complex conjugate pairs, starting from with the number with the lowest real part.
        //If the real part is the same for all the number, order according to the absolute value of the highest imaginary part.
        //Within a pair, the element with negative imaginary part comes first.
        Complex[] Cplxpair(Complex[] complex_vector)
        {

            Complex[] output_complex_vector = new Complex[complex_vector.Length * 2];

            //Order in ascending order
            Complex temp_val;    //Store the real part of the complex number

            //Use Bubble Sort algorithm to order the data in ascending order based on the real-part
            for (int kk = 0; kk < complex_vector.Length; kk++)
            {

                for (int ll = kk + 1; ll < complex_vector.Length; ll++)
                {
                    //If the real parts are different, then sort based on the real part
                    if (complex_vector[kk].Real != complex_vector[ll].Real)
                    {

                        if (complex_vector[kk].Real > complex_vector[ll].Real)
                        {

                            temp_val = complex_vector[kk];
                            complex_vector[kk] = complex_vector[ll];
                            complex_vector[ll] = temp_val;

                        }

                    }

                    //If the real parts are identical, sort based on the imaginary part
                    else
                    {
                        if (complex_vector[kk].Imaginary < complex_vector[ll].Imaginary)
                        {
                            temp_val = complex_vector[kk];
                            complex_vector[kk] = complex_vector[ll];
                            complex_vector[ll] = temp_val;

                        }


                    }

                }

            }

            //Now create the new vector by adding the conjugate. The negative sign is always placed on the top
            int index_complex_vector = 0;
            for (int kk = 0; kk < output_complex_vector.Length; kk += 2)
            {

                output_complex_vector[kk] = complex_vector[index_complex_vector].Real - complex_imag * complex_vector[index_complex_vector].Imaginary;
                output_complex_vector[kk + 1] = complex_vector[index_complex_vector].Real + complex_imag * complex_vector[index_complex_vector].Imaginary;

                index_complex_vector++;

            }

            return output_complex_vector;

        }

        //Step 3: Get N - th order Elliptic analog lowpass prototype
        public void Ellipap(int order_filt, double Rp, double Rs)
        {

            double kc_matlab;
            double kp_matlab;
            double k_matlab;

            double Gp = Math.Pow(10, -Rp / 20);     //passband gain
            double ep = Math.Sqrt(Math.Pow(10, Rp / 10) - 1);    //ripple factors
            double es = Math.Sqrt(Math.Pow(10, Rs / 10) - 1);

            double k1_matlab = ep / es;

            double half_order = Math.Floor((double)order_filt / 2);

            double[] ui_matlab = new double[(int)half_order];    //Initialize the vector with "0" values

            for (int kk = 1; kk < half_order + 1; kk++)
            {

                ui_matlab[kk - 1] = (2 * (double)kk - 1) / order_filt;

            }

            double k_min = Math.Pow(10, -6.0);

            double q_matlab;
            double q1_matlab;

            if (k1_matlab < k_min)
            {

                double[] K_Kprime = Ellipk(k1_matlab);

                q_matlab = Math.Exp(-Math.PI * (K_Kprime[1] / K_Kprime[0]));
                q1_matlab = Math.Pow(q_matlab, (1 / (double)order_filt));

                double temp_k_matlab_I = 0;
                double temp_k_matlab_II = 0;

                //7 is the default number of expansion terms used by Matlab for the Elliptic filter
                for (double kk = 1; kk < 8; kk++)
                {

                    temp_k_matlab_I += Math.Pow(q1_matlab, kk * (kk + 1));
                    temp_k_matlab_II += Math.Pow(q1_matlab, kk * kk);

                }

                k_matlab = 4 * Math.Sqrt(q1_matlab) * Math.Pow(((1 + temp_k_matlab_I) / (1 + 2 * temp_k_matlab_II)), 2);

            }

            else
            {

                kc_matlab = Math.Sqrt(1 - Math.Pow(k1_matlab, 2.0));
                kp_matlab = Math.Pow(kc_matlab, order_filt) * Sne(ui_matlab, kc_matlab);
                k_matlab = Math.Sqrt(1 - Math.Pow(kp_matlab, 2.0));

            }

            int r_matlab = order_filt % 2;

            //Zeros of elliptic rational function
            double[] zi_matlab = Cde(ui_matlab, k_matlab);

            //Filter zeros = poles of elliptic rational function
            Complex[] z_matlab = new Complex[zi_matlab.Length];

            for (int kk = 0; kk < zi_matlab.Length; kk++)
            {

                z_matlab[kk] = complex_imag / (k_matlab * zi_matlab[kk]);

            }

            //Clear vector that is no longer needed
            Array.Clear(zi_matlab, 0, zi_matlab.Length);

            Complex v0_matlab = complex_real - Asne(complex_imag / ep, k1_matlab);

            v0_matlab = v0_matlab * complex_neg_imag / (double)order_filt;

            Complex[] ui_matlab_complex = new Complex[ui_matlab.Length];

            for (int kk = 0; kk < ui_matlab.Length; kk++)
            {

                ui_matlab_complex[kk] = ui_matlab[kk] - complex_imag * v0_matlab;

            }

            Complex[] p_matlab = Cde(ui_matlab_complex, k_matlab);

            for (int kk = 0; kk < p_matlab.Length; kk++)
            {

                p_matlab[kk] *= complex_imag;

            }

            Complex p0_matlab = complex_imag * Sne(complex_imag * v0_matlab, k_matlab);

            double[][] B_matlab = new double[p_matlab.Length][];
            double[][] A_matlab = new double[p_matlab.Length][];

            for (int hh = 0; hh < p_matlab.Length; hh++)
            {

                B_matlab[hh] = new double[3];
                A_matlab[hh] = new double[3];

            }

            for (int kk = 0; kk < 3; kk++)
            {

                for (int ll = 0; ll < B_matlab.Length; ll++)
                {
                    switch (kk)
                    {

                        case 0:

                            B_matlab[ll][kk] = 1;
                            A_matlab[ll][kk] = 1;
                            break;

                        case 1:

                            Complex temp_val = complex_real / z_matlab[ll];
                            B_matlab[ll][kk] = -2 * temp_val.Real;

                            temp_val = complex_real / p_matlab[ll];
                            A_matlab[ll][kk] = -2 * temp_val.Real;
                            break;

                        case 2:

                            B_matlab[ll][kk] = Math.Pow(Complex.Abs(complex_real / z_matlab[ll]), 2);
                            A_matlab[ll][kk] = Math.Pow(Complex.Abs(complex_real / p_matlab[ll]), 2);
                            break;

                    }

                }

            }

            double[][] B_matlab_I = new double[p_matlab.Length + 1][];
            double[][] A_matlab_I = new double[p_matlab.Length + 1][];

            for (int hh = 0; hh < B_matlab_I.Length; hh++)
            {

                B_matlab_I[hh] = new double[3];
                A_matlab_I[hh] = new double[3];

            }

            if (r_matlab == 0)
            {
                for (int kk = 0; kk < p_matlab.Length + 1; kk++)
                {

                    if (kk == 0)
                    {

                        B_matlab_I[kk][0] = Gp;
                        A_matlab_I[kk][0] = 1;

                    }

                    else
                    {

                        for (int ll = 0; ll < 3; ll++)
                        {

                            B_matlab_I[kk][ll] = B_matlab[kk - 1][ll];
                            A_matlab_I[kk][ll] = A_matlab[kk - 1][ll];

                        }

                    }

                }

            }

            else
            {

                for (int kk = 0; kk < p_matlab.Length + 1; kk++)
                {

                    if (kk == 0)
                    {

                        B_matlab_I[kk][0] = 1;

                        A_matlab_I[kk][0] = 1;

                        Complex temp_val = -complex_real / p0_matlab;
                        A_matlab_I[kk][1] = temp_val.Real;

                    }

                    else
                    {

                        for (int ll = 0; ll < 3; ll++)
                        {

                            B_matlab_I[kk][ll] = B_matlab[kk - 1][ll];
                            A_matlab_I[kk][ll] = A_matlab[kk - 1][ll];

                        }

                    }

                }

            }

            Complex[] z_matlab_I = Cplxpair(z_matlab);
            Complex[] p_matlab_I = Cplxpair(p_matlab);

            Complex[] p_matlab_II = new Complex[p_matlab_I.Length + 1];
            if (r_matlab == 1)
            {

                for (int kk = 0; kk < p_matlab_I.Length; kk++)
                {

                    p_matlab_II[kk] = p_matlab_I[kk];

                }

                p_matlab_II[p_matlab_I.Length] = p0_matlab;

                Array.Clear(p_matlab_I, 0, p_matlab_I.Length);
                p_matlab_I = p_matlab_II;

            }

            double H0_matlab = Math.Pow(Gp, 1 - r_matlab);

            Complex prod_p = complex_real;
            Complex prod_z = complex_real;

            for (int kk = 0; kk < p_matlab_I.Length; kk++)
            {

                prod_p *= p_matlab_I[kk];

            }

            for (int kk = 0; kk < z_matlab_I.Length; kk++)
            {

                prod_z *= z_matlab_I[kk];

            }

            k_matlab = Complex.Abs(H0_matlab * prod_p / prod_z);

            //Passing the values of the local variables to the global variables
            z_matlab_ellipap = z_matlab_I;
            p_matlab_ellipap = p_matlab_I;
            k_matlab_ellipap = k_matlab;

        }

        //Intermediate method for Step 4: calculate the coefficients of the polynomial (based on Matlab code)
        public Complex[] Poly(Complex[] temp_array_poly, int col_poly)
        {

            Complex[] coeff_pol_f = new Complex[col_poly + 1];
            coeff_pol_f[0] = 1;

            for (int ll = 0; ll < col_poly; ll++)
            {

                int yy = 0;

                do
                {

                    coeff_pol_f[ll + 1 - yy] = coeff_pol_f[ll + 1 - yy] - temp_array_poly[ll] * coeff_pol_f[ll - yy];
                    yy++;

                } while (yy <= ll);

            }

            return coeff_pol_f;

        }

        //Intermediate method for Step 4: calculate the coefficients of the polynomial (based on Matlab code)
        public double[] Poly_Complex(Complex[] temp_array_poly, int col_poly)
        {

            Complex[] coeff_pol_f_complex = new Complex[col_poly + 1];
            coeff_pol_f_complex[0] = 1;

            double[] coeff_pol_f = new double[col_poly + 1];

            for (int ll = 0; ll < col_poly; ll++)
            {

                int yy = 0;

                do
                {

                    coeff_pol_f_complex[ll + 1 - yy] = coeff_pol_f_complex[ll + 1 - yy] - temp_array_poly[ll] * coeff_pol_f_complex[ll - yy];
                    yy++;

                } while (yy <= ll);

            }

            //Return only the real part
            for (int kk = 0; kk < coeff_pol_f_complex.Length; kk++)
            {

                coeff_pol_f[kk] = coeff_pol_f_complex[kk].Real;

            }

            return coeff_pol_f;

        }

        //Intermediate method for Step 4: calculate the coefficients of the polynomial (based on Matlab code)
        double[] Poly(double[][] temp_array_poly, int col_poly)
        {
            //Extract the eigenvectors and eigenvalues
            Matrix<double> eigval_vec = Matrix<double>.Build.Dense(temp_array_poly.Length, temp_array_poly.Length);

            for (int kk = 0; kk < temp_array_poly.Length; kk++)
            {
                for (int ll = 0; ll < temp_array_poly.Length; ll++)
                {

                    eigval_vec[kk, ll] = temp_array_poly[kk][ll];

                }
            }

            Evd<double> eigen = eigval_vec.Evd();

            //Going through this extra step,because C# was giving me a hard time when 
            //trying to create a vector of complex numbers
            Matrix<Complex> eigval = Matrix<Complex>.Build.Dense(1, eigen.EigenValues.Count);

            for (int kk = 0; kk < eigen.EigenValues.Count; kk++)
            {

                eigval[0, kk] = eigen.EigenValues[kk];

            }


            //Reorganize the eigenvalues in ascending order
            Complex temp_val;
            Complex[] eigval_a = new Complex[eigen.EigenValues.Count];
            for (int kk = 0; kk < eigen.EigenValues.Count; kk++)
            {

                eigval_a[kk] = eigval[0, eigen.EigenValues.Count - 1 - kk];

                if ((kk % 2) == 1 && kk > 0)
                {

                    if (eigval_a[kk - 1].Imaginary < eigval_a[kk].Imaginary)
                    {

                        temp_val = eigval_a[kk - 1];
                        eigval_a[kk - 1] = eigval_a[kk];
                        eigval_a[kk] = temp_val;

                    }

                }

            }

            Complex[] coeff_pol_f_complex = new Complex[col_poly + 1];
            coeff_pol_f_complex[0] = 1;

            double[] coeff_pol_f = new double[col_poly + 1];

            for (int ll = 0; ll < col_poly; ll++)
            {

                int yy = 0;

                do
                {

                    coeff_pol_f_complex[ll + 1 - yy] = coeff_pol_f_complex[ll + 1 - yy] - eigval_a[ll] * coeff_pol_f_complex[ll - yy];
                    yy++;

                } while (yy <= ll);

            }

            //Return only the real part
            for (int kk = 0; kk < coeff_pol_f_complex.Length; kk++)
            {

                coeff_pol_f[kk] = coeff_pol_f_complex[kk].Real;

            }

            return coeff_pol_f;

        }


        //Step 4: Transform to state-space
        public void Zp2ss()
        {

            int order_p = p_matlab_ellipap.Length;
            int order_z = z_matlab_ellipap.Length;

            bool oddpoles_matlab = false;
            bool oddZerosOnly_matlab = false;

            b_matlab = new double[p_matlab_ellipap.Length];
            c_matlab = new double[p_matlab_ellipap.Length];

            a_matlab = new double[p_matlab_ellipap.Length][];

            for (int hh = 0; hh < p_matlab_ellipap.Length; hh++)
            {

                a_matlab[hh] = new double[p_matlab_ellipap.Length];

            }

            d_matlab = 1;

            Complex[] coeff_num;
            Complex[] coeff_den;
            double wn_matlab;

            t_matlab = new double[2][];

            for (int hh = 0; hh < 2; hh++)
            {

                t_matlab[hh] = new double[2];

            }

            //Odd number of poles and zeros
            if (p_matlab_ellipap.Length % 2 == 1 && z_matlab_ellipap.Length % 2 == 1)
            {

                a_matlab[0][0] = p_matlab_ellipap[p_matlab_ellipap.Length - 1].Real;
                b_matlab[0] = 1;
                c_matlab[0] = p_matlab_ellipap[p_matlab_ellipap.Length - 1].Real - z_matlab_ellipap[p_matlab_ellipap.Length - 1].Real;
                d_matlab = 1;

                order_p--;
                order_z--;

                oddpoles_matlab = true;

            }

            //If odd number of poles only
            else if (p_matlab_ellipap.Length % 2 == 1)
            {

                a_matlab[0][0] = p_matlab_ellipap[p_matlab_ellipap.Length - 1].Real;
                b_matlab[0] = 1;
                c_matlab[0] = 1;
                d_matlab = 0;

                order_p--;

                oddpoles_matlab = true;

            }

            //If odd number of zeros only
            else if (z_matlab_ellipap.Length % 2 == 1)
            {

                coeff_num = Poly(z_matlab_ellipap, 2);
                coeff_den = Poly(p_matlab_ellipap, 2);

                wn_matlab = Math.Sqrt(Complex.Abs(p_matlab_ellipap[p_matlab_ellipap.Length - 2] * p_matlab_ellipap[p_matlab_ellipap.Length - 1]));

                if (wn_matlab == 0)
                {

                    wn_matlab = 1;

                }


                t_matlab[0][0] = 1;
                t_matlab[1][1] = 1 / wn_matlab;
                a_matlab[0][0] = t_matlab[0][0] * (-coeff_den[1].Real) / t_matlab[0][0];
                a_matlab[0][1] = t_matlab[0][1] * (-coeff_den[2].Real) / t_matlab[0][1];
                a_matlab[1][0] = t_matlab[1][0] * 1 / t_matlab[1][0];
                a_matlab[1][1] = 0;

                b_matlab[0] = 1 / t_matlab[0][0];

                c_matlab[0] = t_matlab[0][0];
                c_matlab[1] = coeff_num[1].Real;

                oddZerosOnly_matlab = true;

            }

            Complex[] temp_poly_p = new Complex[2];
            Complex[] temp_poly_z = new Complex[2];

            double[][] a1_matlab = new double[2][];
            double[] b1_matlab = new double[2];
            double[] c1_matlab = new double[2];
            double d1_matlab = 1;

            for (int hh = 0; hh < 2; hh++)
            {

                a1_matlab[hh] = new double[2];

            }

            int track_index = 1;
            int j_index;
            while (track_index < order_z)
            {

                for (int rr = track_index - 1; rr < track_index + 1; rr++)
                {

                    temp_poly_p[rr - track_index + 1] = p_matlab_ellipap[rr];
                    temp_poly_z[rr - track_index + 1] = z_matlab_ellipap[rr];

                }

                coeff_num = Poly(temp_poly_z, 2);
                coeff_den = Poly(temp_poly_p, 2);

                for (int kk = 0; kk < coeff_den.Length; kk++)
                {

                    coeff_num[kk] = coeff_num[kk].Real;
                    coeff_den[kk] = coeff_den[kk].Real;

                }

                wn_matlab = Math.Sqrt(Complex.Abs(p_matlab_ellipap[track_index - 1] * p_matlab_ellipap[track_index]));

                if (wn_matlab == 0)
                {

                    wn_matlab = 1;

                }

                t_matlab[0][0] = 1;
                t_matlab[1][1] = 1 / wn_matlab;

                //Since t_matlab is a diagonal matrix, no need to include the values multiplied in the 
                //(1,2) and (2,1) position;
                a1_matlab[0][0] = (t_matlab[0][0] * (-coeff_den[1].Real));
                a1_matlab[0][1] = (t_matlab[1][1] * (-coeff_den[2].Real));
                a1_matlab[1][0] = wn_matlab;
                a1_matlab[1][1] = 0;

                b1_matlab[0] = 1;

                c1_matlab[0] = t_matlab[0][0] * (coeff_num[1].Real - coeff_den[1].Real);
                c1_matlab[1] = t_matlab[1][1] * (coeff_num[2].Real - coeff_den[2].Real);

                d1_matlab = 1;

                if (oddpoles_matlab)
                {

                    j_index = track_index - 1;

                }

                else if (oddZerosOnly_matlab)
                {

                    j_index = track_index;

                }

                else
                {

                    j_index = track_index - 2;

                }

                if (j_index == -1)
                {

                    a_matlab[0][0] = a1_matlab[0][0];
                    a_matlab[0][1] = a1_matlab[0][1];
                    a_matlab[1][0] = a1_matlab[1][0];
                    a_matlab[1][1] = a1_matlab[1][1];

                }

                else
                {

                    //Since b1 has always a zero in the second row, no need to include 
                    //j_index + 2 in the "for" loop
                    for (int kk = 0; kk < j_index + 1; kk++)
                    {

                        a_matlab[j_index + 1][kk] = c_matlab[kk];

                    }

                    a_matlab[j_index + 1][j_index + 1] = a1_matlab[0][0];
                    a_matlab[j_index + 1][j_index + 2] = a1_matlab[0][1];
                    a_matlab[j_index + 2][j_index + 1] = a1_matlab[1][0];
                    a_matlab[j_index + 2][j_index + 2] = a1_matlab[1][1];

                }

                b_matlab[j_index + 1] = b1_matlab[0] * d_matlab;
                b_matlab[j_index + 2] = b1_matlab[1] * d_matlab;

                c_matlab[j_index + 1] = c1_matlab[0];
                c_matlab[j_index + 2] = c1_matlab[1];

                d_matlab *= d1_matlab;

                track_index += 2;

            }

            while (track_index < order_p)
            {

                for (int rr = track_index - 1; rr < track_index + 1; rr++)
                {

                    temp_poly_p[rr - track_index + 1] = p_matlab_ellipap[rr];

                }

                coeff_den = Poly(temp_poly_p, 2);

                for (int kk = 0; kk < coeff_den.Length; kk++)
                {

                    coeff_den[kk] = coeff_den[kk].Real;

                }

                wn_matlab = Math.Sqrt(Complex.Abs(p_matlab_ellipap[track_index - 1] * p_matlab_ellipap[track_index]));

                if (wn_matlab == 1)
                {

                    wn_matlab = 0;

                }

                t_matlab[0][0] = 1;
                t_matlab[1][1] = 1 / wn_matlab;

                //Since t_matlab is a diagonal matrix, no need to include the values multiplied in the 
                //(1,2) and (2,1) position;
                a1_matlab[0][0] = (t_matlab[0][0] * (-coeff_den[1].Real));
                a1_matlab[0][1] = (t_matlab[1][1] * (-coeff_den[2].Real));
                a1_matlab[1][0] = wn_matlab;
                a1_matlab[1][1] = 0;

                b1_matlab[0] = 1;

                c1_matlab[0] = 0;   //t_matlab[0][0];
                c1_matlab[1] = t_matlab[1][1];

                d1_matlab = 0;

                if (oddpoles_matlab)
                {

                    j_index = track_index - 1;

                }

                else if (oddZerosOnly_matlab)
                {

                    j_index = track_index;

                }

                else
                {

                    j_index = track_index - 2;

                }

                if (j_index == -1)
                {

                    a_matlab[0][0] = a1_matlab[0][0];
                    a_matlab[0][1] = a1_matlab[0][1];
                    a_matlab[1][0] = a1_matlab[1][0];
                    a_matlab[1][1] = a1_matlab[1][1];

                    c_matlab[0] = c1_matlab[0];
                    c_matlab[1] = c1_matlab[1];

                }

                else
                {

                    //Since b1 has always a zero in the second row, no need to include 
                    //j_index + 2 in the "for" loop
                    for (int kk = 0; kk < j_index + 1; kk++)
                    {

                        a_matlab[j_index + 1][kk] = c_matlab[kk];
                        c_matlab[kk] = d1_matlab * c_matlab[kk];

                    }

                    a_matlab[j_index + 1][j_index + 1] = a1_matlab[0][0];
                    a_matlab[j_index + 1][j_index + 2] = a1_matlab[0][1];
                    a_matlab[j_index + 2][j_index + 1] = a1_matlab[1][0];
                    a_matlab[j_index + 2][j_index + 2] = a1_matlab[1][1];

                    c_matlab[j_index + 1] = c1_matlab[0];
                    c_matlab[j_index + 2] = c1_matlab[1];

                }

                b_matlab[j_index + 1] = b1_matlab[0] * d_matlab;
                b_matlab[j_index + 2] = b1_matlab[1] * d_matlab;

                d_matlab *= d1_matlab;

                track_index += 2;

            }

            for (int kk = 0; kk < c_matlab.Length; kk++)
            {

                c_matlab[kk] *= k_matlab_ellipap.Real;

            }

            d_matlab *= k_matlab_ellipap.Real;

        }

        //Step 5: Use Bilinear transformation to find discrete equivalent
        public void Bilinear(Matrix<double> a_arma_f, Matrix<double> b_arma_f, Matrix<double> c_arma_f, Matrix<double> d_arma_f, double fs_f, int type_filt_f, int temp_dim_arr_matr)
        {

            double t_arma;
            double r_arma;

            t1_arma = Matrix<double>.Build.Dense(temp_dim_arr_matr, temp_dim_arr_matr);
            t2_arma = Matrix<double>.Build.Dense(temp_dim_arr_matr, temp_dim_arr_matr);
            ad_arma = Matrix<double>.Build.Dense(temp_dim_arr_matr, temp_dim_arr_matr);
            bd_arma = Matrix<double>.Build.Dense(temp_dim_arr_matr, 1);
            cd_arma = Matrix<double>.Build.Dense(1, temp_dim_arr_matr);
            dd_arma = Matrix<double>.Build.Dense(1, 1);

            Matrix<double> t1_arma_identity = Matrix<double>.Build.DenseIdentity(temp_dim_arr_matr, temp_dim_arr_matr);
            Matrix<double> t2_arma_identity = Matrix<double>.Build.DenseIdentity(temp_dim_arr_matr, temp_dim_arr_matr);

            t_arma = (1 / fs_f);
            r_arma = Math.Sqrt(t_arma);
            t1_arma = t1_arma_identity + a_arma_f * t_arma * 0.5; //t1_arma.eye() 
            t2_arma = t2_arma_identity - a_arma_f * t_arma * 0.5;
            ad_arma = t1_arma * t2_arma.PseudoInverse();
            bd_arma = (t_arma / r_arma) * t2_arma.Solve(b_arma_f);
            cd_arma = (r_arma * c_arma_f) * (t2_arma.PseudoInverse());
            dd_arma = (c_arma_f * t2_arma.PseudoInverse()) * b_arma_f * (t_arma / 2) + d_arma_f;

        }

        //Extract the zeros of the state-space system
        public void Sss_zeros()
        {

            Matrix<double> temp_sss = ad_arma - bd_arma * dd_arma.PseudoInverse() * cd_arma;

            Evd<double> eigen = temp_sss.Evd();

            //Going through this extract step,because C# was giving me a hard time when 
            //trying to create a vector of complex numbers
            Complex[] eigval = new Complex[eigen.EigenValues.Count];

            for (int kk = 0; kk < eigen.EigenValues.Count; kk++)
            {

                eigval[kk] = eigen.EigenValues[kk];

            }


            num_filt = Poly_Complex(eigval, eigval.Length);

        }

        //--------------------------------------------------------------------------------------------------//
        //--------------------------------------------------------------------------------------------------//
        //Calculate the coefficients of the band pass filter
        public double[][] Lp2bp(int order_filt, double Rp, double Rs, double W_f1, double W_f2)
        {

            //Clean up the global variables for a new analysis
            if (!(save_filt_coeff == null))
            {

                Array.Clear(z_matlab_ellipap, 0, z_matlab_ellipap.Length);
                Array.Clear(p_matlab_ellipap, 0, p_matlab_ellipap.Length);
                k_matlab_ellipap = 0;

                Array.Clear(save_filt_coeff, 0, save_filt_coeff.Length);

                Array.Clear(num_filt, 0, num_filt.Length);
                Array.Clear(den_filt, 0, den_filt.Length);

            }

            int type_filt = 0;

            //Step 1: get analog, pre - warped frequencies
            Freq_pre_wrapped(type_filt, W_f1, W_f2);

            //Step 2: convert to low-pass prototype estimate
            Wn_f1_Wn_f2(type_filt, u_f1, u_f2);

            //Step 3: Get N - th order Elliptic analog lowpass prototype
            Ellipap(order_filt, Rp, Rs);

            //Step 4: Transform to state-space
            Zp2ss();

            //Copy the values of the matrix/arrays into "arma" matrix/array in order to compute the pseudo-inverse of the matrix and other matrix operations
            Matrix<double> a_arma = Matrix<double>.Build.Dense(2 * a_matlab.Length, 2 * a_matlab.Length);

            Matrix<double> a_arma_p_eye = Matrix<double>.Build.DenseIdentity(a_matlab.Length, a_matlab.Length);
            Matrix<double> a_arma_n_eye = -Matrix<double>.Build.DenseIdentity(a_matlab.Length, a_matlab.Length);

            Matrix<double> b_arma = Matrix<double>.Build.Dense(2 * b_matlab.Length, 1);
            Matrix<double> c_arma = Matrix<double>.Build.Dense(1, 2 * c_matlab.Length);
            Matrix<double> d_arma = Matrix<double>.Build.Dense(1, 1);

            //Tranform from low-pass to band-pass
            double q_matlab = Wn / Bw;

            for (int kk = 0; kk < 2 * a_matlab.Length; kk++)
            {

                if (kk < a_matlab.Length)
                {

                    b_arma[kk, 0] = b_matlab[kk] * Wn / q_matlab;
                    c_arma[0, kk] = c_matlab[kk];

                }

                for (int ll = 0; ll < 2 * a_matlab.Length; ll++)
                {

                    if (kk < a_matlab.Length)
                    {

                        if (ll < a_matlab.Length)

                        {

                            a_arma[kk, ll] = Wn * a_matlab[kk][ll] / q_matlab;

                        }

                        else
                        {

                            a_arma[kk, ll] = Wn * a_arma_p_eye[kk, ll - a_matlab.Length];

                        }
                    }

                    else
                    {
                        if (ll < a_matlab.Length)
                        {

                            a_arma[kk, ll] = Wn * a_arma_n_eye[kk - a_matlab.Length, ll];

                        }
                    }

                }

            }

            d_arma[0, 0] = d_matlab;

            int dim_matrix = 2 * a_matlab.Length;

            //Clean some memory
            Array.Clear(a_matlab, 0, a_matlab.Length);
            Array.Clear(b_matlab, 0, b_matlab.Length);
            Array.Clear(c_matlab, 0, c_matlab.Length);

            //Step 5: Use Bilinear transformation to find discrete equivalent
            Bilinear(a_arma, b_arma, c_arma, d_arma, fs, type_filt, dim_matrix);

            double[][] ad_arma_d = new double[ad_arma.RowCount][];
            for (int kk = 0; kk < ad_arma.RowCount; kk++)
            {

                ad_arma_d[kk] = new double[ad_arma.ColumnCount];

            }

            for (int kk = 0; kk < ad_arma.RowCount; kk++)
            {

                for (int ll = 0; ll < ad_arma.ColumnCount; ll++)
                {

                    ad_arma_d[kk][ll] = ad_arma[kk, ll];

                }

            }

            den_filt = Poly(ad_arma_d, dim_matrix);

            //Step 6: Extract the zeros from the State-Space Model
            Sss_zeros();

            //Multiply the numerator by the gain, which is found to be the first Markov Parameter, which is "dd_arma"
            for (int kk = 0; kk < num_filt.Length; kk++)
            {

                num_filt[kk] *= dd_arma[0, 0];

            }

            //Insert zeros, if necessary at the numerator
            double[] num_filt_zeros = new double[num_filt.Length + den_filt.Length - num_filt.Length];
            if (den_filt.Length - num_filt.Length > 0)
            {

                for (int kk = den_filt.Length - num_filt.Length; kk < num_filt_zeros.Length; kk++)
                {

                    num_filt_zeros[kk] = num_filt[kk - den_filt.Length + num_filt.Length];

                }

                num_filt = num_filt_zeros;

            }

            //Save numerator and denominator into the final matrix
            save_filt_coeff = new double[2][];

            for (int hh = 0; hh < 2; hh++)
            {

                save_filt_coeff[hh] = new double[num_filt.Length];

            }

            for (int kk = 0; kk < num_filt.Length; kk++)
            {

                save_filt_coeff[0][kk] = num_filt[kk];
                save_filt_coeff[1][kk] = den_filt[kk];

            }

            return save_filt_coeff;

        }
        //--------------------------------------------------------------------------------------------------//
        //--------------------------------------------------------------------------------------------------//


        //--------------------------------------------------------------------------------------------------//
        //--------------------------------------------------------------------------------------------------//
        //Calculate the coefficients of the band stop filter
        public double[][] Lp2bs(int order_filt, double Rp, double Rs, double W_f1, double W_f2)
        {

            //Clean up the global variables for a new analysis
            if (!(save_filt_coeff == null))
            {

                Array.Clear(z_matlab_ellipap, 0, z_matlab_ellipap.Length);
                Array.Clear(p_matlab_ellipap, 0, p_matlab_ellipap.Length);
                k_matlab_ellipap = 0;

                Array.Clear(save_filt_coeff, 0, save_filt_coeff.Length);

                Array.Clear(num_filt, 0, num_filt.Length);
                Array.Clear(den_filt, 0, den_filt.Length);

            }

            int type_filt = 1;

            //Step 1: get analog, pre - warped frequencies
            Freq_pre_wrapped(type_filt, W_f1, W_f2);

            //Step 2: convert to low-pass prototype estimate
            Wn_f1_Wn_f2(type_filt, u_f1, u_f2);

            //Step 3: Get N - th order Elliptic analog lowpass prototype
            Ellipap(order_filt, Rp, Rs);

            //Step 4: Transform to state-space
            Zp2ss();

            //Copy the values of the matrix/arrays into "arma" matrix/array in order to compute the pseudo-inverse of the matrix and other matrix operations
            Matrix<double> a_arma = Matrix<double>.Build.Dense(2 * a_matlab.Length, 2 * a_matlab.Length);

            Matrix<double> a_arma_p_eye = Matrix<double>.Build.DenseIdentity(a_matlab.Length, a_matlab.Length);
            Matrix<double> a_arma_n_eye = -Matrix<double>.Build.DenseIdentity(a_matlab.Length, a_matlab.Length);

            Matrix<double> b_arma = Matrix<double>.Build.Dense(2 * b_matlab.Length, 1);
            Matrix<double> c_arma = Matrix<double>.Build.Dense(1, 2 * c_matlab.Length);
            Matrix<double> d_arma = Matrix<double>.Build.Dense(1, 1);

            Matrix<double> a_arma_pinv = Matrix<double>.Build.Dense(2 * a_matlab.Length, 2 * a_matlab.Length);
            Matrix<double> b_arma_temp = Matrix<double>.Build.Dense(2 * a_matlab.Length, 1);
            Matrix<double> c_arma_temp = Matrix<double>.Build.Dense(1, 2 * a_matlab.Length);


            //Tranform from low-pass to band-stop
            double q_matlab = Wn / Bw;

            int dim_matrix = 2 * a_matlab.Length;

            //Copy the values into arma vectors;
            for (int gg = 0; gg < a_matlab.Length; gg++)
            {
                b_arma_temp[gg, 0] = b_matlab[gg];
                c_arma_temp[0, gg] = c_matlab[gg];

                for (int hh = 0; hh < a_matlab.Length; hh++)
                {

                    a_arma[gg, hh] = a_matlab[gg][hh];

                }

            }

            a_arma_pinv = a_arma.PseudoInverse();

            int kk;
            int ll;

            d_arma = d_matlab - c_arma_temp * a_arma_pinv * b_arma_temp;
            c_arma = c_arma_temp * a_arma_pinv;
            b_arma = -(Wn * a_arma_pinv * b_arma_temp) / q_matlab;

            for (kk = 0; kk < a_matlab.Length; kk++)
            {

                for (ll = 0; ll < a_matlab.Length; ll++)
                {

                    a_arma[kk, ll] = Wn * a_arma_pinv[kk, ll] / q_matlab;
                    a_arma[kk, ll + a_matlab.Length] = Wn * a_arma_p_eye[kk, ll];

                    a_arma[kk + a_matlab.Length, ll] = Wn * a_arma_n_eye[kk, ll];

                }

            }

            //Clean some memory
            Array.Clear(a_matlab, 0, a_matlab.Length);
            Array.Clear(b_matlab, 0, b_matlab.Length);
            Array.Clear(c_matlab, 0, c_matlab.Length);

            //Step 5: Use Bilinear transformation to find discrete equivalent
            Bilinear(a_arma, b_arma, c_arma, d_arma, fs, type_filt, dim_matrix);

            double[][] ad_arma_d = new double[ad_arma.RowCount][];
            for (int gg = 0; gg < ad_arma.RowCount; gg++)
            {

                ad_arma_d[gg] = new double[ad_arma.ColumnCount];

            }

            for (int gg = 0; gg < ad_arma.RowCount; gg++)
            {

                for (int hh = 0; hh < ad_arma.ColumnCount; hh++)
                {

                    ad_arma_d[gg][hh] = ad_arma[gg, hh];

                }

            }

            den_filt = Poly(ad_arma_d, dim_matrix);

            //Step 6: Extract the zeros from the State-Space Model
            Sss_zeros();

            //Multiply the numerator by the gain, which is found to be the first Markov Parameter, which is "dd_arma"
            for (int gg = 0; gg < num_filt.Length; gg++)
            {

                num_filt[gg] *= dd_arma[0, 0];

            }

            //Insert zeros, if necessary at the numerator
            double[] num_filt_zeros = new double[num_filt.Length + den_filt.Length - num_filt.Length];
            if (den_filt.Length - num_filt.Length > 0)
            {

                for (int gg = den_filt.Length - num_filt.Length; gg < num_filt_zeros.Length; gg++)
                {

                    num_filt_zeros[gg] = num_filt[gg - den_filt.Length + num_filt.Length];

                }

                num_filt = num_filt_zeros;

            }

            //Save numerator and denominator into the final matrix
            save_filt_coeff = new double[2][];

            for (int hh = 0; hh < 2; hh++)
            {

                save_filt_coeff[hh] = new double[num_filt.Length];

            }

            for (int gg = 0; gg < num_filt.Length; gg++)
            {

                save_filt_coeff[0][gg] = num_filt[gg];
                save_filt_coeff[1][gg] = den_filt[gg];

            }

            return save_filt_coeff;

        }
        //--------------------------------------------------------------------------------------------------//
        //--------------------------------------------------------------------------------------------------//



        //--------------------------------------------------------------------------------------------------//
        //--------------------------------------------------------------------------------------------------//
        //Calculate the coefficients of the high pass filter
        public double[][] Lp2hp(int order_filt, double Rp, double Rs, double W_f2)
        {

            //Clean up the global variables for a new analysis
            if (!(save_filt_coeff == null))
            {

                Array.Clear(z_matlab_ellipap, 0, z_matlab_ellipap.Length);
                Array.Clear(p_matlab_ellipap, 0, p_matlab_ellipap.Length);
                k_matlab_ellipap = 0;

                Array.Clear(save_filt_coeff, 0, save_filt_coeff.Length);

                Array.Clear(num_filt, 0, num_filt.Length);
                Array.Clear(den_filt, 0, den_filt.Length);

            }

            int type_filt = 2;

            //Step 1: get analog, pre - warped frequencies
            Freq_pre_wrapped(type_filt, W_f2, 0);

            //Step 2: convert to low-pass prototype estimate
            Wn_f1_Wn_f2(type_filt, u_f1, u_f2);

            //Step 3: Get N - th order Elliptic analog lowpass prototype
            Ellipap(order_filt, Rp, Rs);

            //Step 4: Transform to state-space
            Zp2ss();

            //Copy the values of the matrix/arrays into "arma" matrix/array in order to compute the pseudo-inverse of the matrix and other matrix operations
            Matrix<double> a_arma = Matrix<double>.Build.Dense(a_matlab.Length, a_matlab.Length);
            Matrix<double> b_arma = Matrix<double>.Build.Dense(b_matlab.Length, 1);
            Matrix<double> c_arma = Matrix<double>.Build.Dense(1, c_matlab.Length);
            Matrix<double> d_arma = Matrix<double>.Build.Dense(1, 1);

            for (int kk = 0; kk < a_matlab.Length; kk++)
            {

                b_arma[kk, 0] = b_matlab[kk];
                c_arma[0, kk] = c_matlab[kk];

                for (int ll = 0; ll < a_matlab.Length; ll++)
                {

                    a_arma[kk, ll] = a_matlab[kk][ll];

                }

            }

            //Tranform from low-pass to high-pass
            d_arma = d_matlab - c_arma * a_arma.PseudoInverse() * b_arma;
            c_arma = c_arma * a_arma.PseudoInverse();
            b_arma = -Wn * (a_arma.PseudoInverse() * b_arma);
            a_arma = Wn * a_arma.PseudoInverse();

            int dim_matrix = a_matlab.Length;

            //Clean some memory
            Array.Clear(a_matlab, 0, a_matlab.Length);
            Array.Clear(b_matlab, 0, b_matlab.Length);
            Array.Clear(c_matlab, 0, c_matlab.Length);

            //Step 5: Use Bilinear transformation to find discrete equivalent
            Bilinear(a_arma, b_arma, c_arma, d_arma, fs, type_filt, dim_matrix);

            double[][] ad_arma_d = new double[ad_arma.RowCount][];
            for (int gg = 0; gg < ad_arma.RowCount; gg++)
            {

                ad_arma_d[gg] = new double[ad_arma.ColumnCount];

            }

            for (int gg = 0; gg < ad_arma.RowCount; gg++)
            {

                for (int hh = 0; hh < ad_arma.ColumnCount; hh++)
                {

                    ad_arma_d[gg][hh] = ad_arma[gg, hh];

                }

            }

            den_filt = Poly(ad_arma_d, dim_matrix);

            //Step 6: Extract the zeros from the State-Space Model
            Sss_zeros();

            //Multiply the numerator by the gain, which is found to be the first Markov Parameter, which is "dd_arma"
            for (int kk = 0; kk < num_filt.Length; kk++)
            {

                num_filt[kk] *= dd_arma[0, 0];

            }

            //Insert zeros, if necessary at the numerator
            double[] num_filt_zeros = new double[num_filt.Length + den_filt.Length - num_filt.Length];
            if (den_filt.Length - num_filt.Length > 0)
            {

                for (int gg = den_filt.Length - num_filt.Length; gg < num_filt_zeros.Length; gg++)
                {

                    num_filt_zeros[gg] = num_filt[gg - den_filt.Length + num_filt.Length];

                }

                num_filt = num_filt_zeros;

            }

            //Save numerator and denominator into the final matrix
            save_filt_coeff = new double[2][];

            for (int hh = 0; hh < 2; hh++)
            {

                save_filt_coeff[hh] = new double[num_filt.Length];

            }

            for (int gg = 0; gg < num_filt.Length; gg++)
            {

                save_filt_coeff[0][gg] = num_filt[gg];
                save_filt_coeff[1][gg] = den_filt[gg];

            }

            return save_filt_coeff;

        }
        //--------------------------------------------------------------------------------------------------//
        //--------------------------------------------------------------------------------------------------//



        //--------------------------------------------------------------------------------------------------//
        //--------------------------------------------------------------------------------------------------//
        //Calculate the coefficients of the low pass filter
        public double[][] Lp2lp(int order_filt, double Rp, double Rs, double W_f1)
        {

            //Clean up the global variables for a new analysis
            if (!(save_filt_coeff == null))
            {

                Array.Clear(z_matlab_ellipap, 0, z_matlab_ellipap.Length);
                Array.Clear(p_matlab_ellipap, 0, p_matlab_ellipap.Length);
                k_matlab_ellipap = 0;

                Array.Clear(save_filt_coeff, 0, save_filt_coeff.Length);

                Array.Clear(num_filt, 0, num_filt.Length);
                Array.Clear(den_filt, 0, den_filt.Length);

            }

            int type_filt = 3;

            //Step 1: get analog, pre - warped frequencies
            Freq_pre_wrapped(type_filt, 0, W_f1);

            //Step 2: convert to low-pass prototype estimate
            Wn_f1_Wn_f2(type_filt, u_f1, u_f2);

            //Step 3: Get N - th order Elliptic analog lowpass prototype
            Ellipap(order_filt, Rp, Rs);

            //Step 4: Transform to state-space
            Zp2ss();

            //Transform lowpass to lowpass (step not included in the zp2ss)
            for (int kk = 0; kk < a_matlab.Length; kk++)
            {

                for (int ll = 0; ll < a_matlab.Length; ll++)
                {

                    a_matlab[kk][ll] *= Wn;

                }

            }

            for (int kk = 0; kk < b_matlab.Length; kk++)
            {

                b_matlab[kk] *= Wn;

            }


            //Copy the values of the matrix/arrays into "arma" matrix/array in order to compute the pseudo-inverse of the matrix and other matrix operations
            Matrix<double> a_arma = Matrix<double>.Build.Dense(a_matlab.Length, a_matlab.Length);
            Matrix<double> b_arma = Matrix<double>.Build.Dense(b_matlab.Length, 1);
            Matrix<double> c_arma = Matrix<double>.Build.Dense(1, c_matlab.Length);
            Matrix<double> d_arma = Matrix<double>.Build.Dense(1, 1);

            for (int kk = 0; kk < a_matlab.Length; kk++)
            {
                b_arma[kk, 0] = b_matlab[kk];
                c_arma[0, kk] = c_matlab[kk];

                for (int ll = 0; ll < a_matlab.Length; ll++)
                {

                    a_arma[kk, ll] = a_matlab[kk][ll];

                }

            }

            d_arma[0, 0] = d_matlab;

            int dim_matrix = a_matlab.Length;

            //Clean some memory
            Array.Clear(a_matlab, 0, a_matlab.Length);
            Array.Clear(b_matlab, 0, b_matlab.Length);
            Array.Clear(c_matlab, 0, c_matlab.Length);

            //Step 5: Use Bilinear transformation to find discrete equivalent
            Bilinear(a_arma, b_arma, c_arma, d_arma, fs, type_filt, dim_matrix);

            double[][] ad_arma_d = new double[ad_arma.RowCount][];
            for (int gg = 0; gg < ad_arma.RowCount; gg++)
            {

                ad_arma_d[gg] = new double[ad_arma.ColumnCount];

            }

            for (int gg = 0; gg < ad_arma.RowCount; gg++)
            {

                for (int hh = 0; hh < ad_arma.ColumnCount; hh++)
                {

                    ad_arma_d[gg][hh] = ad_arma[gg, hh];

                }

            }

            den_filt = Poly(ad_arma_d, dim_matrix);

            //Step 6: Extract the zeros from the State-Space Model
            Sss_zeros();

            //Multiply the numerator by the gain, which is found to be the first Markov Parameter, which is "dd_arma"
            for (int gg = 0; gg < num_filt.Length; gg++)
            {

                num_filt[gg] *= dd_arma[0, 0];

            }

            //Insert zeros, if necessary at the numerator
            double[] num_filt_zeros = new double[num_filt.Length + den_filt.Length - num_filt.Length];
            if (den_filt.Length - num_filt.Length > 0)
            {

                for (int gg = den_filt.Length - num_filt.Length; gg < num_filt_zeros.Length; gg++)
                {

                    num_filt_zeros[gg] = num_filt[gg - den_filt.Length + num_filt.Length];

                }

                num_filt = num_filt_zeros;

            }

            //Save numerator and denominator into the final matrix
            save_filt_coeff = new double[2][];

            for (int hh = 0; hh < 2; hh++)
            {

                save_filt_coeff[hh] = new double[num_filt.Length];

            }

            for (int gg = 0; gg < num_filt.Length; gg++)
            {

                save_filt_coeff[0][gg] = num_filt[gg];
                save_filt_coeff[1][gg] = den_filt[gg];

            }

            return save_filt_coeff;

        }
        //--------------------------------------------------------------------------------------------------//
        //--------------------------------------------------------------------------------------------------//

        //--------------------------------------------------------------------------------------------------//
        //--------------------------------------------------------------------------------------------------//
        //Check the stability of the filter
        public bool Check_stability_iir(double[][] coeff_filt)
        {
            bool stability_flag = true;

            //Coefficients need to be organized in ascending order
            double[] temp_coeff_den = new double[coeff_filt[1].Length];
            for (int kk = 0; kk < coeff_filt[1].Length; kk++)
            {

                temp_coeff_den[kk] = coeff_filt[1][coeff_filt[1].Length - 1 - kk];

            }

            Complex[] roots_den = FindRoots.Polynomial(temp_coeff_den);

            double[] magnitude_roots_den = new double[roots_den.Length];

            for (int kk = 0; kk < roots_den.Length; kk++)
            {

                magnitude_roots_den[kk] = Complex.Abs(roots_den[kk]);

                if (magnitude_roots_den[kk] >= 1)
                {

                    stability_flag = false;
                    break;
                }

            }

            return stability_flag;

        }
        //--------------------------------------------------------------------------------------------------//
        //--------------------------------------------------------------------------------------------------//

        //--------------------------------------------------------------------------------------------------//
        //--------------------------------------------------------------------------------------------------//
        //Filter the data by using the Direct-Form II Transpose, as explained in the Matlab documentation
        //Filter the data by using the Direct-Form II Transpose, as explained in the Matlab documentation
        public double[] Filter_Data(double[][] coeff_filt, double[] pre_filt_signal)
        {

            double[] filt_signal = new double[pre_filt_signal.Length];
            Array.Clear(filt_signal, 0, filt_signal.Length);

            double[][] w_val = new double[coeff_filt[0].Length][];


            for (int ff = 0; ff < coeff_filt[0].Length; ff++)
            {

                w_val[ff] = new double[pre_filt_signal.Length];

            }


            //Convolution product to filter the data
            for (int kk = 0; kk < pre_filt_signal.Length; kk++)
            {

                if (kk == 0)
                {

                    filt_signal[kk] = pre_filt_signal[kk] * coeff_filt[0][0];


                    for (int ww = 1; ww < coeff_filt[0].Length; ww++)
                    {

                        w_val[ww - 1][kk] = pre_filt_signal[kk] * coeff_filt[0][ww] - filt_signal[kk] * coeff_filt[1][ww];


                    }

                }

                else
                {

                    filt_signal[kk] = pre_filt_signal[kk] * coeff_filt[0][0] + w_val[0][kk - 1];

                    for (int ww = 1; ww < coeff_filt[0].Length; ww++)
                    {

                        w_val[ww - 1][kk] = pre_filt_signal[kk] * coeff_filt[0][ww] + w_val[ww][kk - 1] - filt_signal[kk] * coeff_filt[1][ww];

                        if (ww == coeff_filt[0].Length - 1)
                        {

                            w_val[ww - 1][kk] = pre_filt_signal[kk] * coeff_filt[0][ww] - filt_signal[kk] * coeff_filt[1][ww];

                        }

                    }

                }



            }

            return filt_signal;

        }

    }

}
