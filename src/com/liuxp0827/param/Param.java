package com.liuxp0827.param;

import com.liuxp0827.constant.constant;
import com.liuxp0827.math.Jmath;
import com.liuxp0827.utils.JvprException;
import com.liuxp0827.utils.QuoteInteger;
import com.liuxp0827.waveIO.WavInfo;

import java.util.ArrayList;

/**
 * Created by liuxiangpeng on 2016/11/15.
 */
public class Param {
    private FilterBank filterBank;
    private Mfcc mfcc;
    private float[] cepLifterWinSize;
    private double[] hammingWinSize;    // vector of the hamming window
    private int warpWinLength;          // warping window size
    private float[] warpTable;          // warping probability table

    public Param() {
        this.warpWinLength = 300;
    }

    public Mfcc GetMfcc() {
        return this.mfcc;
    }

    // Initialize the filter bank info struct. User should
    //  call this function before calling wav2MFCC().
    // - Arguments -
    //   sampleRate : Sample rate (samples per second)
    //   frameRate : Frame length in ms
    //   filterBankSize : Number of filter banks
    //   lowCutFreq : Low cut-off frequency (in Hz) default = -2
    //   highCutFreq : High cut-off frequency (in Hz)default = -1
    public void InitFBank(int sampleRate, int frameRate, int filterBankSize) throws JvprException {
        this.InitFBank2(sampleRate, frameRate, filterBankSize, -2, -1);
    }

    public void InitFBank2(int sampleRate, int frameRate, int filterBankSize, int lowCutFreq, int highCutFreq) throws JvprException {
        float melLowCutFreq, melHighCutFreq, melStep, curFreq, fa, fb, fc;

        if (lowCutFreq >= highCutFreq || highCutFreq > (sampleRate >> 1) || lowCutFreq > (sampleRate >> 1)) {
            throw new JvprException("Low and High cut off frequencies set incorrectly");
        }

        // check number of filter bancks
        if (this.mfcc != null) {
            if (this.mfcc.mfccOrder + 1 > this.filterBank.filterBankSize) {
                throw new JvprException("param order nb greater than filter bank nb");
            }
        }

        this.filterBank = new FilterBank();

        // given by arguments
        this.filterBank.sampleRate = sampleRate;
        this.filterBank.frameSize = (int) ((float) sampleRate * (float) frameRate * 1e-3);
        this.filterBank.filterBankSize = filterBankSize;

        // calculated from arguments
        this.filterBank.fttSize = 2;
        while (this.filterBank.frameSize > this.filterBank.fttSize) {
            this.filterBank.fttSize <<= 1;
        }

        int fttIndex = this.filterBank.fttSize >> 1;
        this.filterBank.fttResolution = (float) sampleRate / (float) this.filterBank.fttSize;

        // the low and high cut-off indices
        if (lowCutFreq < 0) {
            this.filterBank.start = 0;
            lowCutFreq = 0;
        } else {
            this.filterBank.start = lowCutFreq / (int) (this.filterBank.fttResolution);
        }

        if (highCutFreq < 0) {
            highCutFreq = sampleRate >> 1;
            this.filterBank.end = fttIndex;
        } else {
            this.filterBank.end = (int) ((float) (highCutFreq) / this.filterBank.fttResolution);
            if (this.filterBank.end > fttIndex) {
                this.filterBank.end = fttIndex;
            }
        }

        this.filterBank.centerFreqs = null;
        this.filterBank.lowerFilterBanksIndex = null;
        this.filterBank.lowerFilterBanksWeight = null;
        this.filterBank.fftRealValue = null;
        this.filterBank.fftComplexValue = null;

        // the center frequencies
        this.filterBank.centerFreqs = new float[filterBankSize + 1];

        melLowCutFreq = this.mel((float) lowCutFreq);
        melHighCutFreq = this.mel((float) highCutFreq);
        melStep = (melHighCutFreq - melLowCutFreq) / (float) (filterBankSize + 1);
        this.filterBank.centerFreqs[0] = (float) (lowCutFreq); // the zero index is the low cut-off

        for (int i = 1; i <= filterBankSize; i++) {
            this.filterBank.centerFreqs[i] = melLowCutFreq + melStep * (float) (i);
            this.filterBank.centerFreqs[i] = this.freq(this.filterBank.centerFreqs[i]);
        }

        // lower channel indices
        this.filterBank.lowerFilterBanksIndex = new short[fttIndex];

        int ichan = 0;
        for (int i = 0; i < fttIndex; i++) {
            curFreq = (float) (i) * this.filterBank.fttResolution;
            if (i < this.filterBank.start || i > this.filterBank.end) {
                this.filterBank.lowerFilterBanksIndex[i] = -1;
            } else {
                while (ichan <= filterBankSize && this.filterBank.centerFreqs[ichan] <= curFreq) {
                    ichan++;
                }
                this.filterBank.lowerFilterBanksIndex[i] = (short) (ichan - 1);
            }
        }

        // lower channel weights
        this.filterBank.lowerFilterBanksWeight = new float[fttIndex];
        for (int i = 0; i < fttIndex; i++) {
            curFreq = (float) (i) * this.filterBank.fttResolution;
            if (this.filterBank.lowerFilterBanksIndex[i] == -1) {
                this.filterBank.lowerFilterBanksWeight[i] = 0.0f;
            } else {
                if ((int) (this.filterBank.lowerFilterBanksIndex[i]) < filterBankSize) {
                    fc = this.filterBank.centerFreqs[this.filterBank.lowerFilterBanksIndex[i] + 1];
                } else {
                    fc = (float) (highCutFreq);
                }

                fa = 1 / (this.filterBank.centerFreqs[this.filterBank.lowerFilterBanksIndex[i]] - fc);
                fb = -fa * fc;
                this.filterBank.lowerFilterBanksWeight[i] = fa * curFreq + fb;
            }
        }

        // alloc memory for data buffer
//        this.filterBank.fftRealValue = new [this.filterBank.fttSize];
//        this.filterBank.fftComplexValue = new double[this.filterBank.fttSize];
        this.filterBank.fftRealValue = new ArrayList<Double>();
        this.filterBank.fftComplexValue = new ArrayList<Double>();

        // the defaults
        this.filterBank.isLogFBChannels = true;
        this.filterBank.isUsePower = false;
        this.filterBank.isPreEmphasize = true;
        this.filterBank.isUseHamming = true;
    }


    // Initialize the mfcc info struct. User should call
//  this function before calling wav2MFCC().
// - Arguments -
//      iOrder : MFCC order (except the 0th)
//  fFrameRate : MFCC frame rate in ms

    public void InitMfcc(int iOrder, float fFrmRate) throws JvprException {

        this.mfcc = new Mfcc();

        if (this.filterBank != null) {
            if (iOrder + 1 > this.filterBank.filterBankSize) {
                throw new JvprException("param order nb greater than filter bank nb");
            }
        }

        this.mfcc.mfccOrder = iOrder;  // mfcc order, except the 0th
        this.mfcc.FrameRate = fFrmRate; // frame rate in ms

        // the defaults
        this.mfcc.dynamicWinSize = 2;
        this.mfcc.isFilter = false;
        this.mfcc.IsStatic = true;
        this.mfcc.IsDynamic = true;
        this.mfcc.IsAcce = false;
        this.mfcc.IsLiftCepstral = true;
        this.mfcc.CepstralLifter = 22;

    }

    // wave -> mfcc.
// - Arguments -
//       data : buffer for wave data sequence
//      fParam : buffer for storing the converted parameters,
//               memory is alloced within the function, so the
//               CALLER is RESPONSIBLE to free the memory.
//        col : width of the param vector
//        row : length of the param vector sequence

    public void Wav2Mfcc(float[] data, WavInfo wavinfo, ArrayList<Float> fParam, QuoteInteger col, QuoteInteger row) throws JvprException {

        if (this.mfcc.IsZeroGlobalMean) {
            this.IsZeroGlobalMean(data, wavinfo.Length);
        }

        if (this.mfcc.IsDBNorm) {
            this.dBNorm(data, wavinfo.Length);
        }

        if (this.filterBank == null || this.mfcc == null) {
            throw new JvprException("Filter bank info and MFCC info not initialized");
        }

        QuoteInteger width = new QuoteInteger();
        int iFrameRate;
        int fttIndex = this.filterBank.fttSize >> 1;
        float melfloor = 1.0f;
        float[] fstatic;

        // calculate number of rows (frames)
        iFrameRate = (int) (1e-3 * (float) (this.mfcc.FrameRate) * (float) (this.filterBank.sampleRate));

        if (iFrameRate > this.filterBank.frameSize) {
            throw new JvprException("Sample point equal to zero");
        }

        row.setValue((int) ((wavinfo.Length - (this.filterBank.frameSize - iFrameRate)) / (iFrameRate)));

        // buffer for raw static params (include the 0th coef)
        width.setValue(this.mfcc.mfccOrder + 1);

        fstatic = new float[row.getValue() * width.getValue()];

        // buffer for filter banks

        for (int i = 0; i < row.getValue(); i++) {
//            this.filterBank.fftRealValue = new double[this.filterBank.fttSize];
//            this.filterBank.fftComplexValue = new double[this.filterBank.fttSize];
            this.filterBank.fftRealValue = new ArrayList<Double>(this.filterBank.fttSize);
            this.filterBank.fftComplexValue = new ArrayList<Double>(this.filterBank.fttSize);

            for (int j = 0; j < this.filterBank.frameSize; j++) {
                this.filterBank.fftRealValue.add(j, (double) (data[i * iFrameRate + j]));
            }

            // Do pre-emphasis
            if (this.filterBank.isPreEmphasize) {
                this.preEmphasise(this.filterBank.fftRealValue, this.filterBank.frameSize);
            }

            // Do hamming
            if (this.filterBank.isUseHamming) {
                this.doHamming(this.filterBank.fftRealValue, this.filterBank.frameSize);
            }

            // take fft
            Jmath.FFT(this.filterBank.fftRealValue, this.filterBank.fftComplexValue, this.filterBank.fttSize);

            double[] filterBank = new double[this.filterBank.filterBankSize];
            for (int j = 0; j < fttIndex; j++) {
                this.filterBank.fftRealValue.set(j
                        , this.filterBank.fftRealValue.get(j) * this.filterBank.fftRealValue.get(j)
                                + this.filterBank.fftComplexValue.get(j) * this.filterBank.fftComplexValue.get(j));
//                this.filterBank.fftRealValue[j] =
//                        this.filterBank.fftRealValue[j] * this.filterBank.fftRealValue[j] +
//                                this.filterBank.fftComplexValue[j] * this.filterBank.fftComplexValue[j];
            }

            //	Differential Power Spectrum
            if (this.filterBank.isUsePower && this.mfcc.IsDiffPowerSpectrum) {
                this.DPSCC(fttIndex);
            }

            // use power or amp
            if (!this.filterBank.isUsePower) {
                for (int j = 0; j < fttIndex; j++) {
                    this.filterBank.fftRealValue.set(j, Math.sqrt(this.filterBank.fftRealValue.get(j)));
                }
            }

            // Predictive Differential Amplitude Spectrum
            if (!this.filterBank.isUsePower && this.mfcc.IsPredDiffAmpSpetrum) {
                this.PDASCC(fttIndex);
            }

            // accumulate filter banks
            for (int j = 0; j < fttIndex; j++) {
                if (this.filterBank.lowerFilterBanksIndex[j] < 0) {
                    continue;
                }

                // accumulate the lower bank
                if (this.filterBank.lowerFilterBanksIndex[j] != 0) {
                    filterBank[this.filterBank.lowerFilterBanksIndex[j] - 1] +=
                            this.filterBank.fftRealValue.get(j) * (double) (this.filterBank.lowerFilterBanksWeight[j]);
                }
                // accumulate the upper bank
                if ((this.filterBank.lowerFilterBanksIndex[j]) < this.filterBank.filterBankSize) {
                    filterBank[this.filterBank.lowerFilterBanksIndex[j]] +=
                            this.filterBank.fftRealValue.get(j) * (double) (1 - this.filterBank.lowerFilterBanksWeight[j]);
                }
            }

            // take logs
            if (this.filterBank.isLogFBChannels) {
                for (int j = 0; j < this.filterBank.filterBankSize; j++) {
                    if (filterBank[j] >= (melfloor)) {
                        filterBank[j] = Math.log(filterBank[j]);
                    } else {
                        filterBank[j] = Math.log((double) melfloor);
                    }
                }
            }

            // take dct
            if (!this.mfcc.isFilter) {
                filterBank = Jmath.DCT(filterBank, width);


                // Liftering
                if (this.mfcc.IsLiftCepstral) {
                    this.liftCepstral(filterBank);
                }
            }

            // copy data
            for (int j = 0; j < width.getValue(); j++) {
                fstatic[i * width.getValue() + j] = (float) (filterBank[j]);
            }
        }

        if (this.mfcc.IsFeatWarping) {
            this.warping(fstatic, width.getValue(), row, width.getValue(), (int) this.mfcc.FeatWarpWinSize);

        }

        if (this.mfcc.IsRasta) {
            this.rastaFiltering(fstatic, width.getValue(), row, width.getValue());
        }

        if (this.mfcc.IsEnergyNorm) {
            this.energyNorm(fstatic, width.getValue(), row.getValue());
        }

        this.static2Full(fstatic, width, row, fParam);

        col.setValue(width.getValue());
    }

    // Calculate the requested parameters from the raw static coefficients.
//  This function convert raw static coef. to conjunct coef. (0th, static,
//  delta, acce or any conbinations of them) specified in pmfccinfo.
//  In calculating deltas, some heading and tailing frames may be discarded
//  in this function. The caller is responsible to to free the memory of
//  the returned dbuffer.
// - Arguments -
//     fStatic : buffer for raw static coef. (include the 0th)
//        col : width of the raw static coef., and is changed
//               by this function to be the width of the conjunct
//               params.
//        row : frames of the raw static coef. and is recalculated
//               in this function to be the actual number of frames
//                of the conjunct params.
    public void static2Full(float[] fstatic, QuoteInteger col, QuoteInteger row, ArrayList<Float> retFParam) throws JvprException {
        int width = col.getValue();
        int iSOff, iDOff, ipt = 0;
        float[] fdelta = null, facce = null;
        Float[] fParam;

        if (this.mfcc == null) {
            retFParam.clear();
            retFParam = null;
            throw new JvprException("MFCC info not initialized");
        }

        // take deltas from statics
        if (this.mfcc.IsDynamic) {
            fdelta = new float[width * row.getValue()];
            this.doDelta(fdelta, fstatic, row, width);


            iSOff = this.mfcc.dynamicWinSize;
            iDOff = 0;
        }

        // take accelerations from deltaes
        if (this.mfcc.IsAcce) {
            facce = new float[width * row.getValue()];
            this.doDelta(facce, fdelta, row, width);


            iSOff = 2 * this.mfcc.dynamicWinSize;
            iDOff = this.mfcc.dynamicWinSize;
        }

        iSOff = 0;
        iDOff = 0;

        // calculate the actual width of the conjunct parameter
        col.setValue(0);

        if (this.mfcc.IsStatic) {
            col.add(this.mfcc.mfccOrder);
        }

        if (this.mfcc.IsDynamic) {
            col.add(this.mfcc.mfccOrder);
        }

        if (this.mfcc.IsAcce) {
            col.add(this.mfcc.mfccOrder);
        }

        // prepare for parameter buffer
        fParam = new Float[col.getValue() * row.getValue()];

        for (int i = 0; i < row.getValue(); i++) {

            if (this.mfcc.IsStatic) {
                for (int j = 1; j < width; j++) {
                    fParam[ipt] = fstatic[(i + iSOff) * width + j];
                    ipt++;
                }
            }

            if (this.mfcc.IsDynamic) {
                for (int j = 1; j < width; j++) {
                    fParam[ipt] = fdelta[(i + iDOff) * width + j];
                    ipt++;
                }
            }

            if (this.mfcc.IsAcce) {
                for (int j = 1; j < width; j++) {
                    fParam[ipt] = facce[i * width + j];
                    ipt++;
                }
            }
        }

        for (int i = 0; i < fParam.length; i++) {
            retFParam.set(i, fParam[i]);
        }
    }

    // ------------- Cepstral Mean Substraction & Variance Normalisation ------------
    // This function normalizes the mfcc feature parameters into a Guassian
    // distribute,which can reduce the influence of channel.
    //  fParam      : buffer which stored feature parameters
    //	iVecsize    : size of a feature vector which stored parameter
    //  iVecNum     : number of feature vectors
    public void FeatureNorm(float[][] fParam, int iVecSize, int iVecNum) throws JvprException {
        if (iVecSize <= 0) {
            throw new JvprException("Dimension of GMM less than zero");
        }

        if (iVecNum <= 0) {
            throw new JvprException("Nb of frames less than zero");
        }

        float[] cmsMean = new float[iVecSize];
        float[] cmsStdv = new float[iVecSize];
        float tempMean = 0f, tempStdv = 0f;

        //Get the average value of the mV
        for (int i = 0; i < iVecSize / 2; i++) {
            for (int j = 0; j < iVecNum; j++) {
                tempMean += fParam[j][i];
                tempStdv += fParam[j][i] * fParam[j][i];
            }

            cmsMean[i] = tempMean / (float) iVecNum;

            //Get the standard deviations
            cmsStdv[i] = tempStdv / (float) iVecNum;
            cmsStdv[i] -= cmsMean[i] * cmsMean[i];

            if (cmsStdv[i] <= 0) {
                cmsStdv[i] = 1.0f;
            } else {
                cmsStdv[i] = (float) Math.sqrt(cmsStdv[i]);
            }

            tempMean = 0;
            tempStdv = 0;
        }

        //subtract the average value
        for (int i = 0; i < iVecSize / 2; i++) {
            for (int j = 0; j < iVecNum; j++) {
                fParam[j][i] = (fParam[j][i] - cmsMean[i]) / cmsStdv[i];
            }
        }

    }

    public void FeatureNorm2(float[] fParam, int iVecSize, int iVecNum) throws JvprException {
        if (iVecSize <= 0) {
            throw new JvprException("Dimension of GMM less than zero");
        }

        if (iVecNum <= 0) {
            throw new JvprException("Nb of frames less than zero");
        }

        float[] cmsMean = new float[iVecSize];
        float[] cmsStdv = new float[iVecSize];
        float tempMean = 0f, tempStdv = 0f;

        //Get the average value of the mV
        for (int i = 0; i < iVecSize / 2; i++) {
            for (int j = 0; j < iVecNum; j++) {
                tempMean += fParam[j * iVecSize + i];
                tempStdv += fParam[j * iVecSize + i] * fParam[j * iVecSize + i];
            }

            cmsMean[i] = tempMean / (float) (iVecNum);

            //Get the standard deviations
            cmsStdv[i] = tempStdv / (float) (iVecNum);
            cmsStdv[i] -= cmsMean[i] * cmsMean[i];

            if (cmsStdv[i] <= 0) {
                cmsStdv[i] = 1.0f;
            } else {
                cmsStdv[i] = (float) (Math.sqrt(cmsStdv[i]));
            }

            tempMean = 0;
            tempStdv = 0;
        }

        //subtract the average value
        for (int i = 0; i < iVecSize / 2; i++) {
            for (int j = 0; j < iVecNum; j++) {
                fParam[j * iVecSize + i] = (fParam[j * iVecSize + i] - cmsMean[i]) / cmsStdv[i];
            }
        }
    }

    //------------- Decibel Normalization -------------------------------------------
    public void dBNorm(float[] sampleBuffer, long sampleCount) {
        float sampleMax = -Float.MAX_VALUE;
        for (int i = 0; i < sampleCount; i++) {
            if (sampleBuffer[i] > sampleMax) {
                sampleMax = sampleBuffer[i];
            }
        }

        for (int i = 0; i < sampleCount; i++) {
            sampleBuffer[i] = (float) (Math.pow(10, constant.DB / 20.0)) * (float) (Math.pow(2, 15) - 1) * sampleBuffer[i] / sampleMax;
        }
    }

    //------------- IsZeroGlobalMean --------------------------------------------------
    public void IsZeroGlobalMean(float[] data, long sampleCount) {
        float mean = 0.0f;
        for (int i = 0; i < sampleCount; i++) {
            mean += data[i];
        }
        mean /= (float) sampleCount;

        for (int i = 0; i < sampleCount; i++) {
            float y = data[i] - mean;
            if (y > 32767) {
                y = 32767;
            }
            if (y < -32767) {
                y = -32767;
            }
            if (y > 0) {
                data[i] = (float) ((short) (y + 0.5));
            } else {
                data[i] = (float) ((short) (y - 0.5));
            }
        }
    }

    //------------- Norm static feature ---------------------------------------------
    private void energyNorm(float[] p_FeatBuf, int p_nVecSize, int p_nFrameNum) {
        float maxe, mine;
        float[] ft = p_FeatBuf;
        int index = 0;
        maxe = ft[index];

        for (int i = 0; i < p_nFrameNum; i++) {
            if (ft[index] > maxe) {
                maxe = ft[index];
            }

            index += p_nVecSize;
        }

        mine = (maxe - (float) (this.mfcc.SilFloor) * (float) (Math.log(10.0))) / 10.0f;

        for (int i = 0; i < p_nFrameNum; i++) {
            if (ft[index] < mine) {
                mine = ft[index];
            }

            ft[index] = 1.0f - (maxe - ft[index]) * (float) (this.mfcc.EnergyScale);
            p_FeatBuf[index] = 1.0f - (maxe - p_FeatBuf[index]) * (float) (this.mfcc.EnergyScale);
            index += p_nVecSize;
        }
    }

    //------------- Differential Power Spectrum -------------------------------------
    public void DPSCC(int pointNB) {
        int fttIndex = pointNB;
        for (int j = 0; j < fttIndex; j++) {
            if (j < fttIndex - 1) {
                this.filterBank.fftRealValue.set(j, Math.abs(this.filterBank.fftRealValue.get(j) - this.filterBank.fftRealValue.get(j + 1)));
            } else {
                this.filterBank.fftRealValue.set(j, 0.0);
            }
        }
    }

    //------------- Predictive Differential Amplitude Spectrum ----------------------
    public void PDASCC(int pointNB) {
        int fttIndex = pointNB;

        //	1.预测差分
        int WINLEN = 6;
        double[] damplitude = new double[fttIndex];
        for (int j = 0; j < fttIndex; j++) {
            double dmax = -Double.MAX_VALUE;
            for (int w = 0; j + w < fttIndex && w <= WINLEN; w++) {
                double dsin = Math.sin((w * constant.PI) / 2 * WINLEN);
                double dcur = this.filterBank.fftRealValue.get(j + w) * dsin;
                if (dcur > dmax) {
                    dmax = dcur;
                }
            }
            damplitude[j] = dmax;
        }

        double alpha = 1.05;
        double[] dDright = new double[fttIndex];
        double[] dDleft = new double[fttIndex];
        for (int j = 0; j < fttIndex - 1; j++) {
            if (damplitude[j] > this.filterBank.fftRealValue.get(j) && damplitude[j + 1] < this.filterBank.fftRealValue.get(j + 1)) {
                dDright[j] = this.filterBank.fftRealValue.get(j) - alpha * this.filterBank.fftRealValue.get(j + 1);
            } else if (damplitude[j] <= this.filterBank.fftRealValue.get(j) && damplitude[j + 1] >= this.filterBank.fftRealValue.get(j + 1)) {
                dDright[j] = alpha * this.filterBank.fftRealValue.get(j) - this.filterBank.fftRealValue.get(j + 1);
            } else {
                dDright[j] = this.filterBank.fftRealValue.get(j) - this.filterBank.fftRealValue.get(j + 1);
            }
        }

        dDright[fttIndex - 1] = 0.0;

        for (int j = fttIndex - 1; j > 0; j--) {
            if (damplitude[j] < this.filterBank.fftRealValue.get(j) && damplitude[j - 1] < this.filterBank.fftRealValue.get(j - 1)) {
                dDleft[j] = this.filterBank.fftRealValue.get(j) - alpha * this.filterBank.fftRealValue.get(j - 1);
            } else if (damplitude[j] >= this.filterBank.fftRealValue.get(j) && damplitude[j - 1] >= this.filterBank.fftRealValue.get(j - 1)) {
                dDleft[j] = alpha * this.filterBank.fftRealValue.get(j) - this.filterBank.fftRealValue.get(j - 1);
            } else {
                dDleft[j] = this.filterBank.fftRealValue.get(j) - this.filterBank.fftRealValue.get(j - 1);
            }
        }

        if (damplitude[0] < this.filterBank.fftRealValue.get(0)) {
            dDleft[0] = (1.0 - alpha) * this.filterBank.fftRealValue.get(0);
        } else {
            dDleft[0] = (alpha - 1.0) * this.filterBank.fftRealValue.get(0);
        }

        // 2.累积过程

        double[] left = new double[fttIndex];
        double[] right = new double[fttIndex];
        for (int i = 1; i < fttIndex; i++) {
            left[i] = left[i - 1] + dDleft[i - 1];
        }

        for (int i = fttIndex - 2; i >= 0; i--) {
            right[i] = right[i + 1] + dDright[i + 1];
        }

        for (int i = 0; i < fttIndex; i++) {
            this.filterBank.fftRealValue.set(i, (left[i] + right[i]) / 2.0);
        }
    }

    //------------- Feature Warping -------------------------------------------------
// nWinSize=300
// vSize=nStep=特征维数
// nInNum 输入特征帧数
// data 输入特征
// nOutNum 输出特征帧数
    private void warping(float[] data, int vSize, QuoteInteger nInNum, int nStep, int nWinSize) throws JvprException {

        this.createWarpTable();
        int nOutNum = nInNum.getValue() - nWinSize;
        if (nOutNum <= 0) {
            throw new JvprException("nOutNum can not <= 0");
        }

        float[] warpBuf = new float[nOutNum * nStep];
        int warpFrmNo = 0;
        float[] pDataIvt = new float[nOutNum * nStep];

        for (int i = 0; i < nStep; i++) {
            //		var dst []float32 = make([]float32, nOutNum*nStep+i*nInNum, nOutNum*nStep+i*nInNum)
            for (int j = 0; j < nInNum.getValue(); j++) {
                //			dst[j] = data[j*nStep+i]
                if (j < nOutNum * nStep) {
                    pDataIvt[j] = data[j * nStep + i];
                }
            }
        }

        int halfwin = nWinSize >> 1;
        float[] minus_res = new float[nWinSize];

        for (int i = halfwin; i + halfwin < nInNum.getValue(); i++) {
            for (int k = 0; k < vSize; k++) {

                float[] p = new float[pDataIvt.length + nOutNum * nStep + k * nInNum.getValue() - nOutNum * nStep];
                for (int c = 0; c < pDataIvt.length; c++) {
                    p[c] = pDataIvt[c];
                }


                float curValue = p[i];
                int t = halfwin - i;
                int nIndex = 2 * halfwin - 1;

                for (int m = i - halfwin; m < i; m++) {
                    minus_res[t + m] = p[m] - curValue;
                }

                t = halfwin - i - 1;

                for (int m = i + 1; m < i + halfwin; m++) {
                    minus_res[t + m] = p[m] - curValue;
                }

                int[] ui = new int[nWinSize];
                for (int ii = 0; ii < nWinSize; ii++) {
                    ui[ii] = (int) (minus_res[ii]) & 0x00000000ffffffff;
                }

                for (int mm = 0; mm < 2 * halfwin - 1; mm++) {
                    nIndex -= (ui[mm] >> 31);
                }

                warpBuf[warpFrmNo * nStep + k] = this.warpTable[nIndex];
            }

            for (int k = vSize; k < nStep; k++) {
                warpBuf[warpFrmNo * nStep + k] = data[i * nStep + k];
            }
            warpFrmNo++;
        }

        nInNum.setValue(warpFrmNo);

        System.arraycopy(warpBuf, 0, data, 0, warpFrmNo * nStep);

//        copy(data[:warpFrmNo * nStep],warpBuf[:warpFrmNo * nStep])
    }

//------------- Rasta-filtering -------------------------------------------------
/************************************************************************/
/*
    data : static mfcc
	vSize: order of mfcc
	nNum : frame number
	nStep: order of mfcc
*/

    /************************************************************************/
    private void rastaFiltering(float[] data, int vSize, QuoteInteger nNum, int nStep) throws JvprException {
        if (nNum.getValue() <= 4) {
            throw new JvprException("Not eoungh features for Rasta filtering");
        }

        float[] RastaBuf = new float[nNum.getValue() * vSize];
        for (int i = 0; i < nNum.getValue() - 4; i++) {
            if (i == 0) {
                for (int j = 0; j < vSize; j++) {
                    RastaBuf[i * vSize + j] = 0.1f * (2.0f * data[(i + 4) * nStep + j] + data[(i + 3) * nStep + j] -
                            data[(i + 1) * nStep + j] - 2.0f * data[i * nStep + j]);
                }
            } else {
                for (int j = 0; j < vSize; j++) {
                    RastaBuf[i * vSize + j] = 0.1f * (2.0f * data[(i + 4) * nStep + j] + data[(i + 3) * nStep + j] -
                            data[(i + 1) * nStep + j] - 2.0f * data[i * nStep + j]) +
                            (float) this.mfcc.RastaCoff * RastaBuf[(i - 1) * vSize + j];
                }
            }
        }

        for (int i = 0; i < nNum.getValue() - 4; i++) {
            for (int j = 0; j < vSize; j++) {
                data[i * nStep + j] = RastaBuf[i * vSize + j];
            }
        }

        nNum.sub(4);
    }

    //------------- Create warping talbe --------------------------------------------
    private void createWarpTable() {
        double TableBegin = -10.0, TableEnd = 10.0, presice = 1.0e-5;
        this.warpTable = new float[this.warpWinLength];
        double[] rankBuf = new double[this.warpWinLength];

        for (int i = 0; i < this.warpWinLength; i++) {
            rankBuf[i] = ((double) this.warpWinLength - 0.5 - i) / (double) (this.warpWinLength);
        }

        double integral = 0.0;
        int Index = this.warpWinLength - 1;

        for (double x = TableBegin; x <= TableEnd; x += presice) {
            integral += Math.exp(-x * x / 2.0) / Math.sqrt(2 * constant.PI) * presice;
            if (integral >= rankBuf[Index]) {
                this.warpTable[Index] = (float) x;
                Index--;
                if (Index < 0) {
                    break;
                }
            }
        }
        return;
    }

    // mel -> frequency
    private float freq(float mel) {
        return (float) (700 * (Math.exp((double) (mel) / (double) (1127)) - 1));
    }

    // frequency -> mel
    private float mel(float freq) {
        return (float) (1127 * Math.log(1 + (double) (freq) / (double) (700)));
    }

    // Hamming window
    // - Arguments -
    //  dVector : vector to be windowed
    //     iLen : length of the vector
    private void doHamming(ArrayList<Double> dVector, int iLen) throws JvprException {
        if (dVector == null || dVector.size() < iLen) {
            throw new JvprException("dVector is null or size of dVector < iLen");
        }

        double a;
        if (this.hammingWinSize != null && iLen != this.hammingWinSize.length) {
            this.hammingWinSize = new double[iLen];
        }

        if (this.hammingWinSize == null) {
            this.hammingWinSize = new double[iLen];
            a = (double) 2 * constant.PI / (double) (iLen - 1);
            for (int i = 0; i < iLen; i++) {
                this.hammingWinSize[i] = 0.54 - 0.46 * Math.cos(a * (double) (i));
            }
        }

        for (int i = 0; i < iLen; i++) {
//            dVector[i] *= this.hammingWinSize[i];
            dVector.set(i, dVector.get(i) * this.hammingWinSize[i]);
        }
    }

    // Pre-emphasize the input signal.
    // - Arguments -
    //           s : pointer to input vector
    //        iLen : length of the input vector
    //        preE : option for the emphasis filter default = 0.97
    private void preEmphasise(ArrayList<Double> s, int iLen) throws JvprException {
        this.preEmphasise2(s, iLen, 0.97);
    }

    private void preEmphasise2(ArrayList<Double> s, int iLen, double preE) throws JvprException {
        if (s.size() < iLen) {
            throw new JvprException("ArrayList s'size %d < iLen %d", s.size(), iLen);
        }

        for (int i = iLen; i >= 2; i--) {
            s.set(i - 1, s.get(i - 1) - s.get(i - 2) * preE);
//            s[i - 1] -= s[i - 2] * preE;
        }
//        s[0] *= 1.0 - preE;
        s.set(0, s.get(0) * (1.0 - preE));
    }

    //  Take delta from the relatively static signal. the deltaes
    //  are appended to the static signal. Therefore the given
    //  buffer must be large enough. Note that some frames will
    //  be discarded in this procedure, the memory is freed at
    //  the same time.
    private void doDelta(float[] fdest, float[] fsource, QuoteInteger iLen, int width) throws JvprException {
        int winSize = this.mfcc.dynamicWinSize;
        if (iLen.getValue() < 2 * winSize + 1) {
            throw new JvprException("iLen = %d less than %d", iLen, 2 * winSize + 1);
        }

        float fnorm = 0.0f, fsum;
        float[] fpback, fpforw;

        if (!this.mfcc.IsPolishDiff) {

            for (int k = 1; k <= winSize; k++) {
                fnorm += (float) (k * k);
            }
        } else {
            for (int k = 1; k <= winSize; k++) {
                fnorm += winSize - k + 1;
            }
        }

        fnorm *= 2;

        for (int i = 0; i < iLen.getValue(); i++) {
            for (int d = 0; d < width; d++) {
                fsum = 0;
                for (int k = 1; k <= winSize; k++) {

                    fpback = new float[fsource.length - d + __max(i - k, 0) * width];
                    System.arraycopy(fsource, d + __max(i - k, 0) * width, fpback, 0, fsource.length - d + __max(i - k, 0) * width);
                    fpforw = new float[fsource.length - d + __min(i + k, iLen.getValue() - 1) * width];
                    System.arraycopy(fsource, d + __min(i + k, iLen.getValue() - 1) * width, fpforw, 0, fsource.length - d + __min(i + k, iLen.getValue() - 1) * width);
                    //  fpback = fsource[d + __max(i - k, 0) * width:]
                    //  fpforw = fsource[d + __min(i + k, iLen.getValue() - 1) * width:]
                    float im;
                    if (!this.mfcc.IsPolishDiff) {
                        im = k;
                    } else {
                        im = (float) (winSize - k + 1) / (float) (k);
                    }
                    fsum = fsum + im * (fpforw[0] - fpback[0]);
                }
                fdest[i * width + d] = fsum / fnorm;
            }
        }
    }

    // Lift the cepstral to the same amplitudes. It should be
    //  called just after the dct procedure and before the
    //  deletion of the 0th coefficient.
    // - Arguments -
    //     fVector : vector to be lifted, length is specified in pmfccinfo.
    private void liftCepstral(double[] dVector) throws JvprException {
        float L;
        if (this.mfcc == null) {
            throw new JvprException("CParam mfcc can not be nil");
        }
        int iLen = this.mfcc.mfccOrder + 1;
        L = this.mfcc.CepstralLifter;

        if (this.cepLifterWinSize == null) {
            this.cepLifterWinSize = new float[iLen];
            for (int i = 0; i < iLen; i++) {
                this.cepLifterWinSize[i] = 1.0f + L / 2.0f * (float) (Math.sin(constant.PI * i / (double) L));
            }
        }

        for (int i = 0; i < iLen; i++) {
            dVector[i] *= (double) (this.cepLifterWinSize[i]);
        }

    }

    private int __max(int x, int y) {
        if (x > y) {
            return x;
        }
        return y;
    }

    private int __min(int x, int y) {
        if (x > y) {
            return y;
        }
        return x;
    }


}
