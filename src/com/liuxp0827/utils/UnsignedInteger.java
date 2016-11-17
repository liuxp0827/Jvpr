package com.liuxp0827.utils;

/**
 * Created by liuxiangpeng on 2016/11/16.
 */
public class UnsignedInteger {


//    public static short getUint8(short s) {
//        return s & 0x00ff;
//    }


    public static int getUint16(int i) {
        return i & 0x0000ffff;
    }


    public static long getUint32(long l) {
        return l & 0x00000000ffffffff;
    }
}
