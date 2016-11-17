package com.liuxp0827.test;

import com.liuxp0827.utils.QuoteInteger;

import java.util.ArrayList;

/**
 * Created by liuxiangpeng on 2016/11/16.
 */
public class Test {
    private double[] hammingWinSize;

    public static void main(String[] args) {
        Test test = new Test();
//        ArrayList<Float> outer = new ArrayList<Float>();
//        test.TestValueCopy(outer);
//        if (outer == null) {
//            System.err.println("outer is null");
//            return;
//        }
//
//        for (int i = 0; i < outer.size(); i++) {
//            System.out.println("outer[" + i + "] : " + outer.get(i));
//        }

//        QuoteInteger i2 = new QuoteInteger(2);
//        test.TestIngeterCopy(i2);
//        System.out.println(i2.getValue());
//        System.out.println(test.hammingWinSize.length);
//        test.TestIngeterCopy(i2);
//        System.out.println(test.hammingWinSize.length);
//        test.TestIngeterCopy(i2);
//        System.out.println(test.hammingWinSize == null);

        ArrayList<Double> da = new ArrayList<Double>(100);
        da.add(0, 1.0);
        da.add(1, 2.0);
        da.add(2, 3.0);
        da.add(3, 4.0);
        da.add(4, 5.0);
        System.out.println(da.get(0));
        da.set(0, 1.1);

        System.out.println(da.size());
        System.out.println(da.get(0));
    }


    public void TestValueCopy(ArrayList<Float> fArray) {
        float[] inner = new float[10];
        for (int i = 0; i < inner.length; i++) {
            inner[i] = (i + 1) * 1.0f;
            fArray.add(i, (i + 1) * 1.0f);
        }

        for (int i = 0; i < fArray.size(); i++) {
            System.out.println("fArray[" + i + "] : " + fArray.get(i));
        }

        fArray.clear();
        fArray = null;
    }

    public void TestIngeterCopy(QuoteInteger ii) {
        ii.setValue(111);
        if (this.hammingWinSize == null) {
            this.hammingWinSize = new double[100];
        } else if (this.hammingWinSize != null && this.hammingWinSize.length != 101) {
            this.hammingWinSize = new double[101];
        } else {
            this.hammingWinSize = null;
        }
    }
}


