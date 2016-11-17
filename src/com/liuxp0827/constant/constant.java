package com.liuxp0827.constant;

/**
 * Created by liuxiangpeng on 2016/11/15.
 */
public class constant {
    // Speaker Recognition
    public static final int REL_FACTOR = 16;        // 自适应因子
    public static final int MAXLOP = 1;             // 自适应次数
    public static final int BIT_PER_SAMPLE = 16;
    public static final int SAMPLERATE = 16000;     // 采样率
    public static final int MIN_FRAMES = 300;
    public static final float DB = -3.0f;           // 归一化分贝量
    public static final double LOGZERO = (-1.0E10); // ~log(0)
    public static final double LSMALL = (-0.5E10);  // log values < LSMALL are set to LOGZERO

    public static final double PI = 3.14159265358979;

    public static final double VAR_FLOOR = 0.005;   // Variance floor, make sure the variance to be large enough	(old:0.005)
    public static final double VAR_CEILING = 10.0;  // Variance ceiling
    public static final int MAX_LOOP = 10;

    // Wave
    public static final int VOC_BLOCK_LEN = 1600;
    public static final int MIN_VOC_ENG = 500;
    public static final int SHRT_MAX = 35535;


    public static final double DLOG2PAI = 1.837877066;


    //	Feature  -------------------------
    public static final int LOW_CUT_OFF = 250;
    public static final int HIGH_CUT_OFF = 3800;
    public static final int FILTER_BANK_SIZE = 24;
    public static final int FRAME_LENGTH = 20;
    public static final int FRAME_SHIFTt = 10;
    public static final int MFCC_ORDER = 16;
    public static final boolean BSTATIC = true;
    public static final boolean BDYNAMIC = true;
    public static final boolean BACCE = false;

    public static final boolean CMSVN = true;
    public static final boolean DBNORM = true;
    public static final boolean ZEROGLOBALMEAN = true;
    public static final boolean FEATWARP = false;
    public static final boolean DIFPOL = false;
    public static final boolean DPSCC = false;
    public static final boolean PDASCC = false;
    public static final boolean ENERGYNORM = false;
    public static final boolean RASTA = false;

    public static final int SIL_FLOOR = 50;
    public static final int ENERGY_SCALE = 19;
    public static final int FEATURE_WARPING_WIN_SIZE = 300;
    public static final float RASTA_COFF = 0.94f;

}
