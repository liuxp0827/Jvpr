package com.liuxp0827.param;

import java.util.ArrayList;

/**
 * Created by liuxiangpeng on 2016/11/15.
 */
public class FilterBank {

    protected int sampleRate;         // sample rate (samples / second)
    protected int frameSize;          // frame size (in samples)
    protected int fttSize;            // fft size (2^N)
    protected float fttResolution;    // fft resolution ( (SampRate/2) / (fftN/2) = SampRate/fftN )
    protected int filterBankSize;     // number of filterbanks
    protected int start;
    protected int end;                // lopass to hipass cut-off fft indices

    protected float[] centerFreqs;            // array[1..filterBankSize+1] of centre freqs
    protected short[] lowerFilterBanksIndex;  // array[1..fttSize/2] of lower fbank index
    protected float[] lowerFilterBanksWeight; // array[1..fttSize/2] of lower fbank weighting
    protected ArrayList<Double> fftRealValue;          // array[1..fttSize] of fft bins (real part)
    protected ArrayList<Double> fftComplexValue;       // array[1..fttSize] of fft bins (image part)
//    protected double[] fftRealValue;          // array[1..fttSize] of fft bins (real part)
//    protected double[] fftComplexValue;       // array[1..fttSize] of fft bins (image part)

    protected boolean isUsePower;         // use power rather than magnitude (d: false)
    protected boolean isLogFBChannels;    // log filterbank channels (d: true)
    protected boolean isPreEmphasize;     // pre emphasize (d: true)
    protected boolean isUseHamming;       // Use Hamming window rather than rectangle window (d: true)
}
