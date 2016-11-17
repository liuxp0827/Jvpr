package com.liuxp0827.waveIO;

import java.util.ArrayList;

/**
 * Created by liuxiangpeng on 2016/11/16.
 */
public class FmtChunk {
    byte[] fmt;     // format sub-chunk identification  (4 bytes)
    long flength;   // length of format sub-chunk (4 byte integer)
    short format;   // format specifier (2 byte integer)
    short chans;    // number of channels (2 byte integer)
    long sampsRate; // sample rate in Hz (4 byte integer)
    long bpsec;     // bytes per second (4 byte integer)
    short bpsample; // bytes per sample (2 byte integer)
    short bpchan;   // bits per channel (2 byte integer)

    public ArrayList<Byte> Bytes() {
        ArrayList<Byte> buf = new ArrayList<>();
        // fmt
        for (int i = 0; i < this.fmt.length; i++) {
            buf.add(this.fmt[i]);
        }
        // flength
        buf.add((byte) (this.flength & 0xff));
        buf.add((byte) ((this.flength >> 8) & 0xff));
        buf.add((byte) ((this.flength >> 16) & 0xff));
        buf.add((byte) ((this.flength >> 24) & 0xff));
        // format
        buf.add((byte) (this.format & 0xff));
        buf.add((byte) ((this.format >> 8) & 0xff));
        // chans
        buf.add((byte) (this.chans & 0xff));
        buf.add((byte) ((this.chans >> 8) & 0xff));
        // sampsRate
        buf.add((byte) (this.sampsRate & 0xff));
        buf.add((byte) ((this.sampsRate >> 8) & 0xff));
        buf.add((byte) ((this.sampsRate >> 16) & 0xff));
        buf.add((byte) ((this.sampsRate >> 24) & 0xff));
        // bpsec
        buf.add((byte) (this.bpsec & 0xff));
        buf.add((byte) ((this.bpsec >> 8) & 0xff));
        buf.add((byte) ((this.bpsec >> 16) & 0xff));
        buf.add((byte) ((this.bpsec >> 24) & 0xff));
        // bpsample
        buf.add((byte) (this.bpsample & 0xff));
        buf.add((byte) ((this.bpsample >> 8) & 0xff));
        // bpchan
        buf.add((byte) (this.bpchan & 0xff));
        buf.add((byte) ((this.bpchan >> 8) & 0xff));

        return buf;
    }
}
