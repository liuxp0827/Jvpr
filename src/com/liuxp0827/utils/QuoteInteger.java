package com.liuxp0827.utils;

/**
 * Created by liuxiangpeng on 2016/11/16.
 */
public class QuoteInteger {
    int value;

    public QuoteInteger() {

    }

    public QuoteInteger(int value) {
        this.value = value;
    }

    public int getValue() {
        return value;
    }

    public void setValue(int value) {
        this.value = value;
    }

    public void add(int value) {
        this.value += value;
    }

    public void sub(int value) {
        this.value -= value;
    }
}