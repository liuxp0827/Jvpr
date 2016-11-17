package com.liuxp0827.waveIO;

import java.util.ArrayList;

/**
 * Created by liuxiangpeng on 2016/11/16.
 */
public class DataChunk {
    byte[] data; // data sub-chunk identification  (4 bytes)
    long dlength; // length of data sub-chunk (4 byte integer)

    public ArrayList<Byte> Bytes() {
        ArrayList<Byte> buf = new ArrayList<>();

        for (int i = 0; i < this.data.length; i++) {
            buf.add(this.data[i]);
        }

        buf.add((byte) (this.dlength & 0xff));
        buf.add((byte) ((this.dlength >> 8) & 0xff));
        buf.add((byte) ((this.dlength >> 16) & 0xff));
        buf.add((byte) ((this.dlength >> 24) & 0xff));

        return buf;
    }
}
