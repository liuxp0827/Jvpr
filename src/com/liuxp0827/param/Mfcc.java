package com.liuxp0827.param;

/**
 * Created by liuxiangpeng on 2016/11/15.
 */
public class Mfcc {
    protected int mfccOrder;      // MFCC order (except the 0th)
    protected int dynamicWinSize; // length dynamic window when take delta params, default is 2
    protected boolean isFilter;   // output fbank other than MFCC

    public boolean IsStatic;                 // static coefficients or not
    public boolean IsDynamic;                // dynamic coefficients or not
    public boolean IsAcce;                   // acceleration coefficients or not
    public boolean IsLiftCepstral;           // lift the cepstral or not
    public float FrameRate;                  // frame rate in ms
    public float CepstralLifter;             // cepstral lifter. It's invalid when bFBank is set to true
    public boolean IsPolishDiff;             // polish differential formula
    public boolean IsDBNorm;                 // decibel normalization
    public boolean IsDiffPowerSpectrum;      // Differential Power Spectrum
    public boolean IsPredDiffAmpSpetrum;     // Predictive Differential Amplitude Spetrum
    public boolean IsZeroGlobalMean;
    public boolean IsEnergyNorm;
    public short SilFloor;
    public short EnergyScale;
    public boolean IsFeatWarping;
    public short FeatWarpWinSize;
    public boolean IsRasta;
    public double RastaCoff;

    public Mfcc() {

    }
    
}
