package com.liuxp0827.gmm;

import com.liuxp0827.constant.constant;
import com.liuxp0827.utils.JvprException;

import java.io.*;

/**
 * Created by ponpon on 2016/11/16.
 */
public class Gmm {
    public int Frames;                  // number of total frames
    public float[][] FeatureData;       // feature buffer

    public int VectorSize;              // Vector size
    public int Mixtures;                // Mixtures of the GMM
    public double[] deterCovariance;    // determinant of the covariance matrix [mixture]
    public double[] MixtureWeight;      // weight of each mixture[mixture]
    public double[][] Mean;             // mean vector [mixture,dimension]
    public double[][] Covar;            // covariance (diagonal) [mixture,dimension]

    public Gmm() {

    }

    public void DupModel(Gmm gmm) {
        this.Mixtures = gmm.Mixtures;
        this.VectorSize = gmm.VectorSize;
        this.deterCovariance = new double[this.Mixtures];
        this.MixtureWeight = new double[this.Mixtures];
        this.Mean = new double[this.Mixtures][];
        this.Covar = new double[this.Mixtures][];

        for (int i = 0; i < this.Mixtures; i++) {
            this.deterCovariance[i] = gmm.deterCovariance[i];
            this.MixtureWeight[i] = gmm.MixtureWeight[i];
            this.Mean[i] = new double[this.VectorSize];
            this.Covar[i] = new double[this.VectorSize];
            for (int j = 0; j < this.VectorSize; j++) {
                this.Mean[i][j] = gmm.Mean[i][j];
                this.Covar[i][j] = gmm.Covar[i][j];
            }
        }
    }

    public void LoadModel(String filename) throws IOException {
        FileInputStream fin = new FileInputStream(filename);
        DataInputStream reader = new DataInputStream(fin);

        this.Mixtures = reader.readInt();

        this.VectorSize = reader.readInt();

        this.deterCovariance = new double[this.Mixtures];
        this.MixtureWeight = new double[this.Mixtures];
        this.Mean = new double[this.Mixtures][];
        this.Covar = new double[this.Mixtures][];
        for (int i = 0; i < this.Mixtures; i++) {
            this.Mean[i] = new double[this.VectorSize];
            this.Covar[i] = new double[this.VectorSize];
        }

        for (int i = 0; i < this.Mixtures; i++) {
            this.deterCovariance[i] = 0.0;
        }

        for (int i = 0; i < this.Mixtures; i++) {
            this.MixtureWeight[i] = reader.readDouble();
        }

        for (int i = 0; i < this.Mixtures; i++) {
            reader.readDouble(); // not used
            reader.readDouble(); // not used
            reader.readByte(); // not used

            for (int j = 0; j < this.VectorSize; j++) {
                this.Covar[i][j] = reader.readDouble();
                this.deterCovariance[i] += Math.log(this.Covar[i][j]);
            }

            for (int j = 0; j < this.VectorSize; j++) {
                this.Mean[i][j] = reader.readDouble();
            }
        }
        reader.close();
        fin.close();
    }

    public void SaveModel(String filename) throws IOException {
        FileOutputStream fout = new FileOutputStream(filename);
        DataOutputStream writer = new DataOutputStream(fout);

        writer.writeInt(this.Mixtures);
        writer.writeInt(this.VectorSize);

        for (int i = 0; i < this.Mixtures; i++) {
            writer.writeDouble(this.MixtureWeight[i]);

        }

        for (int i = 0; i < this.Mixtures; i++) {
            writer.writeDouble(0.0);    // not used
            writer.writeDouble(0.0);    // not used
            writer.writeByte(0);        // not used

            for (int j = 0; j < this.VectorSize; j++) {
                writer.writeDouble(this.Covar[i][j]);
            }

            for (int j = 0; j < this.VectorSize; j++) {
                writer.writeDouble(this.Mean[i][j]);
            }
        }

        writer.flush();
        writer.close();
        fout.flush();
        fout.close();
    }

    public void CopyFeatureData(Gmm gmm) {

        this.Frames = gmm.Frames;
        this.VectorSize = gmm.VectorSize;
        this.FeatureData = new float[gmm.Frames][];
        for (int i = 0; i < this.Frames; i++) {
            this.FeatureData[i] = new float[this.VectorSize];
        }

        for (int i = 0; i < this.Frames; i++) {
            for (int j = 0; j < this.VectorSize; j++) {
                this.FeatureData[i][j] = gmm.FeatureData[i][j];
            }
        }
    }


    public int EM(int mixtures) throws JvprException {
        double dlogfrmprob, rubbish, lastrubbish;
        double[] dsumgama, dlogmixw, dgama, dmixw;
        double threshold = 1e-5;
        double[][] mean, covar;
        int loop = 0;

        mean = new double[mixtures][];
        covar = new double[mixtures][];
        for (int i = 0; i < mixtures; i++) {
            mean[i] = new double[this.VectorSize];
            covar[i] = new double[this.VectorSize];
        }

        dmixw = new double[mixtures];
        dlogmixw = new double[mixtures];
        dgama = new double[mixtures];
        dsumgama = new double[mixtures];
        rubbish = .0;

        do {

            lastrubbish = rubbish;
            rubbish = .0;
            for (int i = 0; i < mixtures; i++) {

                // speed up
                if (this.MixtureWeight[i] <= 0) {
                    dlogmixw[i] = constant.LOGZERO;
                } else {
                    dlogmixw[i] = Math.log(this.MixtureWeight[i]);
                }

                // clean up temporary values
                dmixw[i] = .0;
                dsumgama[i] = .0;
                for (int j = 0; j < this.VectorSize; j++) {
                    mean[i][j] = .0;
                    covar[i][j] = .0;
                }
            }

//            for (int i = 0; i < this.Frames; i++) {
//                dlogfrmprob = constant.LOGZERO;
//                for (int j = 0; j < mixtures; j++) {
//                    dgama[j] = this.LMixProb(this.FeatureData[i], j);
//                    dgama[j] += dlogmixw[j];
//                    dlogfrmprob = this.LogAdd(dgama[j], dlogfrmprob);
//                }
//
//                rubbish += dlogfrmprob;
//
//                for (int j = 0; j < mixtures; j++) {
//                    dgama[j] -= dlogfrmprob;
//                    dgama[j] = Math.exp(dgama[j]);
//                    dsumgama[j] += dgama[j];
//
//                    // update weights
//                    dmixw[j] += dgama[j];
//                    for (int k = 0; k < this.VectorSize; k++) {
//                        mean[j][k] += dgama[j] * this.FeatureData[i][k];
//                        covar[j][k] += dgama[j] * (double) (this.FeatureData[i][k]) * (double) (this.FeatureData[i][k]);
//                    }
//                }
//            }

            rubbish /= (double) (this.Frames); // rubbish = LLR

            for (int i = 0; i < mixtures; i++) {
                if (dsumgama[i] == .0) {
                    throw new JvprException("error EM");
                }

                this.MixtureWeight[i] = dmixw[i] / (double) (this.Frames);

                for (int j = 0; j < this.VectorSize; j++) {
                    this.Mean[i][j] = mean[i][j] / dsumgama[i];
                    this.Covar[i][j] = covar[i][j] / dsumgama[i];
                    this.Covar[i][j] -= this.Mean[i][j] * this.Mean[i][j];

                    if (this.Covar[i][j] < constant.VAR_FLOOR) {
                        this.Covar[i][j] = constant.VAR_FLOOR;
                    }

                    if (this.Covar[i][j] > constant.VAR_CEILING) {
                        this.Covar[i][j] = constant.VAR_CEILING;
                    }
                }
            }

            for (int i = 0; i < mixtures; i++) {
                this.deterCovariance[i] = .0;
                for (int j = 0; j < this.VectorSize; j++) {
                    this.deterCovariance[i] += Math.log(this.Covar[i][j]);
                }
            }
            loop++;

        } while (loop < constant.MAX_LOOP && Math.abs((rubbish - lastrubbish) / (lastrubbish + 0.01)) > threshold);


        return loop;
    }

}
