package com.liuxp0827.math;

import com.liuxp0827.constant.constant;
import com.liuxp0827.utils.JvprException;
import com.liuxp0827.utils.QuoteInteger;
import com.liuxp0827.utils.UnsignedInteger;

import java.util.ArrayList;

/**
 * Created by liuxiangpeng on 2016/11/16.
 */
public class Jmath {
    private static double[] dct1;
    private static double[] dct2;

    private static double[] fft1;
    private static double[] fft2;

    //------------------- Fast Fourier Transformation ------------------
    // routine of fft
    // - Arguments -
    //      ar  	: pointer to the sequence of the real part
    //      ai 		: pointer to the sequence of the imaginary part
    //      length  : length of the vector, it must be 2^N, otherwise
    //               the vector must be padded with zero.
    public static void FFT(ArrayList<Double> ar, ArrayList<Double> ai, int length) throws JvprException {

        if (ar == null || ar.size() == 0 || ai == null || ai.size() == 0 || length <= 0) {
            throw new JvprException("invalid param");
        }

        int m = (int) (Math.log1p((double) length));
        if (Math.pow(2.0, (double) (m)) != (double) (length)) {
            throw new JvprException("invalid length");
        }

        double temr, temi, x, y, temr1, temi1;
        if (fft1 == null) {
            fft1 = new double[length];
        } else if (fft1.length != length) {
            fft1 = new double[length];
        }

        if (fft2 == null) {
            fft2 = new double[length];

        } else if (fft2.length != length) {
            fft2 = new double[length];
        }

        // fill buffers with precalculated values
        temr = 2 * constant.PI / (double) (length);
        for (int i = 0; i < length; i++) {
            fft1[i] = Math.sin(temr * (double) (i));
            fft2[i] = Math.cos(temr * (double) (i));
        }

        fftSort(ar, ai, length);

        int b = 1;

        for (int l = 1; l <= m; l++) {
            for (int j = 0; j < b; j++) {

                int p = (1 << UnsignedInteger.getUint32(m - l)) * j;
                x = fft2[p];
                y = -fft1[p];
                p = 1 << UnsignedInteger.getUint32(l);

                for (int k = j; k < length; k += p) {
                    temr = ar.get(k + b) * x - ai.get(k + b) * y;
                    temi = ar.get(k + b) * y + ai.get(k + b) * x;
                    temr1 = ar.get(k) + temr;
                    temi1 = ai.get(k) + temi;
                    ar.set(k + b, ar.get(k) - temr);
                    ai.set(k + b, ai.get(k) - temi);
                    ar.set(k, temr1);
                    ai.set(k, temi1);
                }
            }
            b <<= 1;
        }
    }

    // Discrete Cosine Transform
    // - Arguments -
    //		 data : pointer to the real vector
    //      width : only the first "width" elements will
    //                be computed, if width <= 0, all elements
    //                will be computed and returned
    public static double[] DCT(double[] data, QuoteInteger width) {
        int ioffset;
        int length = data.length;
        int length2 = length << 1;
        int[] icfft = new int[]{1, -1};
        double dfactor = Math.sqrt(2.0 / (double) (length));
        double dc0 = 1.0 / Math.sqrt(2.0);
        if (width.getValue() <= 0) {
            width.setValue(length);
        }

        if (dct1 == null) {
            dct1 = new double[length];
        } else if (fft1.length != length) {
            dct1 = new double[length];
        }

        if (dct2 == null) {
            dct2 = new double[length2];

        } else if (fft2.length != length) {
            dct2 = new double[length2];
        }

        // Create Cosine Table & Copy data to the source buff
        for (int i = 0; i < length2; i++) {
            dct2[i] = Math.cos((double) (i) * constant.PI / (double) (length2));
            if (i < length) {
                dct1[i] = data[i];
            }
        }

        //-------------------------
        // Main Procedure
        for (int k = 0; k < width.getValue(); k++) {
            data[k] = 0.0;
            for (int n = 0; n < length; n++) {
                ioffset = ((n << 1) + 1) * k;
                data[k] += dct1[n] * (double) (icfft[(ioffset / length2) % 2]) * dct2[ioffset % length2];
            }
            if (k == 0) {
                data[k] *= dc0;
            }
            data[k] *= dfactor;
        }

        return data;
    }

    // sort routine for FFT
    private static void fftSort(ArrayList<Double> realx, ArrayList<Double> imagex, int length) {


        int j = length >> 1;
        double tmpRealx, tmpImagex;
        for (int i = 1; i < length - 1; i++) {
            if (i < j) {
                {
                    tmpRealx = realx.get(i);
                    realx.set(i, realx.get(j));
                    realx.set(j, tmpRealx);
                }

                {
                    tmpImagex = imagex.get(i);
                    imagex.set(i, imagex.get(j));
                    imagex.set(j, tmpImagex);
                }
            }

            int k = length >> 1;
            while (j >= k) {
                j = j - k;
                k >>= 1;
            }
            j = j + k;
        }
    }
}
