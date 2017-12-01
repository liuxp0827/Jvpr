package com.liuxp0827.waveIO;

import java.util.ArrayList;

/**
 * Created by liuxiangpeng on 2016/11/16.
 */
public class WaveIO {
    WaveChunk waveChunk; // 12 Bytes
    FmtChunk fmtChunk;   // 24 Bytes
    DataChunk dataChunk; // 8 Bytes

    public byte[] Bytes() {

        ArrayList<Byte> buf = new ArrayList<>();

        buf.addAll(this.waveChunk.Bytes());
        buf.addAll(this.fmtChunk.Bytes());
        buf.addAll(this.dataChunk.Bytes());

        byte[] buffer = new byte[buf.size()];

        for (int i = 0; i < buf.size(); i++) {
            buffer[i] = buf.get(i);
        }

        return buffer;
    }

    public WaveIO(){
        this.waveChunk = new WaveChunk();
        this.fmtChunk = new FmtChunk();
        this.dataChunk = new DataChunk();
    }
}
