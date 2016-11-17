package com.liuxp0827.waveIO;

import java.util.ArrayList;

/**
 * Created by liuxiangpeng on 2016/11/16.
 */
public class WaveChunk {
    byte[] riff; // RIFF file identification (4 bytes)
    long length; // length field (4 bytes)
    byte[] wave; // WAVE chunk identification (4 bytes)

    public ArrayList<Byte> Bytes() {
        ArrayList<Byte> buf = new ArrayList<>();

        for (int i = 0; i < this.riff.length; i++) {
            buf.add(this.riff[i]);
        }

        buf.add((byte) (this.length & 0xff));
        buf.add((byte) ((this.length >> 8) & 0xff));
        buf.add((byte) ((this.length >> 16) & 0xff));
        buf.add((byte) ((this.length >> 24) & 0xff));
        for (int i = 0; i < this.wave.length; i++) {
            buf.add(this.wave[i]);
        }
        return buf;
    }
}
