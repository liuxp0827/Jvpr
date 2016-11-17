package com.liuxp0827.feature;

/**
 * Created by liuxiangpeng on 2016/11/15.
 */
public class Feature {
    private int lowCutOff;                      // low cut-off
    private int highCutOff;                     // high cut-off
    private int filterBankSize;                 // # num of filter-bank
    private int frameLength;                    // # frame length
    private int frameShift;                     // 10 # frame shift
    private int mfccOrder;                      // 16 # Mfcc order
    private boolean isStatic;                   // t	# static Mfcc
    private boolean isDynamic;                  // t	# dynamic Mfcc
    private boolean isAcce;                     // f	# acce Mfcc
    private boolean cmsvn;                      // t	# cmsvn
    private boolean isZeroGlobalMean;           // t # zero global mean
    private boolean isDBNorm;                   // t # decibel normalization
    private boolean isDiffPolish;               // f	# polish differential formula
    private boolean isDiffPowerSpectrum;        // f	# differentail power spectrum
    private boolean isPredDiffAmplSpectrum;     // f	# predictive differential amplitude spectrum
    private boolean isEnergyNorm;
    private short silFloor;
    private short energyscale;
    private boolean isFeatWarping;
    private short featWarpWinSize;
    private boolean isRasta;
    private double rastaCoff;

}
