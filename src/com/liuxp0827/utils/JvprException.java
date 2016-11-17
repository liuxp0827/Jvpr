package com.liuxp0827.utils;

/**
 * Created by liuxiangpeng on 2016/11/16.
 */
public class JvprException extends Exception {

    public JvprException(String str) {
        super(str);
    }

    public JvprException(String format, Object... args) {
        super(String.format(format, args));
    }
}
