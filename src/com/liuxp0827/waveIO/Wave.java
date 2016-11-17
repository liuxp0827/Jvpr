package com.liuxp0827.waveIO;

import com.liuxp0827.utils.JvprException;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;

/**
 * Created by liuxiangpeng on 2016/11/16.
 */
public class Wave {

    public static void WaveSave(String filePath, short[] wavData, long sampsRate, long length) throws IOException {
        WaveIO waveIO = new WaveIO();
        File file = new File(filePath);
        FileOutputStream fos = new FileOutputStream(file);


        // waveChunk
        waveIO.waveChunk.riff = new byte[]{
                'R', 'I', 'F', 'F'
        };
        waveIO.waveChunk.wave = new byte[]{
                'W', 'A', 'V', 'E'
        };
        waveIO.waveChunk.length = length * 2 + 36;

        // fmtChunk
        waveIO.fmtChunk.bpchan = 16;
        waveIO.fmtChunk.bpsample = 2;
        waveIO.fmtChunk.bpsec = 32000;
        waveIO.fmtChunk.chans = 1;
        waveIO.fmtChunk.flength = 16;
        waveIO.fmtChunk.fmt = new byte[]{
                'f', 'm', 't', ' '
        };
        waveIO.fmtChunk.format = 1;
        waveIO.fmtChunk.sampsRate = sampsRate;

        // dataChunk
        waveIO.dataChunk.data = new byte[]{
                'd', 'a', 't', 'a'
        };
        waveIO.dataChunk.dlength = length * 2;

        fos.write(waveIO.Bytes());


        int lenOfWav16 = wavData.length;
        byte[] data = new byte[lenOfWav16 * 2];
        for (int i = 0; i < lenOfWav16 * 2 - 1; ) {
            data[i] = (byte) (wavData[i / 2] & 0xff);
            data[i + 1] = (byte) ((wavData[i / 2] >> 8) & 0xff);
            i += 2;
        }

        fos.write(data);
        fos.flush();
        fos.close();
    }

    public short[] WaveLoad(String filePath) throws IOException, JvprException {
        byte[] header = new byte[44];
        int iTotalReaded = 0, iBytesReaded = 0;
        int lengthOfData;
        byte[] cBuff = new byte[0x4000];

        ArrayList<Byte> wavData = new ArrayList();


        FileInputStream fis = new FileInputStream(filePath);
        fis.read(header);

        if (header[0] != 'R' || header[1] != 'I' || header[2] != 'F' || header[3] != 'F') {
            throw new JvprException("wav format error");
        }

        if (header[22] != 1 || header[23] != 0) {
            throw new JvprException("this wave channel is not 1");
        }

        lengthOfData = (int) (header[40]);
        lengthOfData |= (int) (header[41]) << 8;
        lengthOfData |= (int) (header[42]) << 16;
        lengthOfData |= (int) (header[43]) << 24;

        if (lengthOfData <= 0) {
            throw new JvprException("length of wave data is 0");
        }

        while (true) {
            iBytesReaded = fis.read(cBuff);
            if (iTotalReaded >= lengthOfData) {
                break;
            }
            iTotalReaded += iBytesReaded;
            if (iTotalReaded >= lengthOfData) {
                iBytesReaded = iBytesReaded - (iTotalReaded - lengthOfData);
            }

            for (int i = 0; i < iBytesReaded; i++) {
                wavData.add(cBuff[i]);
            }
        }

        fis.close();

        short[] buf = new short[wavData.size() / 2];
        for (int i = 0, j = 0; i < wavData.size() && j < buf.length; ) {
            buf[j] = (short) (wavData.get(i));
            buf[j] |= (short) (wavData.get(i + 1) << 8);
            i += 2;
            j++;
        }
        return buf;
    }
}
