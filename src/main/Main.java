package main;

public class Main {

    public static Double l = Math.PI;
    public static Double T = 10.0;
    
    public static Integer n;
    public static Integer m;
    public static Double h;
    public static Double tau;
    public static Double TH;

    public static Double A;
    public static Double B;
    public static Double C;

    public static Double A_kn;
    public static Double B_kn;
    public static Double C_kn;
    public static Double c_kn;

    public static Double[] alpha;
    public static Double[] betta;
    public static Double[][] u;

    public static Integer nStep = 5;
    public static Integer mStep = 40;

    public static void main(String[] args) {
        System.out.println("=====M1=====");
        for (int i = 5; i <= 35; i += 10) {
            for (int j = 40; j <= 320; j *=2) {
                setParams(i, j);
                execM1();
                execE("m1");
            }
            System.out.println();
        }
        System.out.println("=====M2=====");
        for (int i = 5; i <= 35; i += 10) {
            for (int j = 40; j <= 320; j *=2) {
                setParams(i, j);
                execM2();
                execE("m2");
            }
            System.out.println();
        }
        System.out.println("=====M3=====");
        for (int i = 5; i <= 35; i += 10) {
            for (int j = 40; j <= 320; j *=2) {
                setParams(i, j);
                execM3();
                execE("m3");
            }
            System.out.println();
        }
    }

    public static void setParams(int i, int j) {
        n = i;
        m = j;
        h = l / n;
        tau = T / m;
        TH = tau / (h * h);
        A = TH;
        B = TH;
        C = -(2 * TH + 1);
        A_kn = TH / 2;
        B_kn = TH / 2;
        C_kn = -(TH + 1);
        c_kn = -(1 - TH);
    }

    private static void prepareData() {
        alpha = new Double[n + 1];
        betta = new Double[n + 1];

        u = new Double[n + 1][m + 1];
        for (int i = 0; i <= n; i++) {
            u[i][0] = Math.sin(getXi(i));
        }
        for (int j = 0; j <= m; j++) {
            u[0][j] = getMu(j);
            u[n][j] = u[0][j];
        }
    }

    public static void execM1() {
        prepareData();
        // явная схема
        for (int j = 0; j < m; j++) {
            for (int i = 1; i < n; i++) {
                u[i][j + 1] = TH * (u[i + 1][j] - 2 * u[i][j] + u[i - 1][j]) + getTauGij(i, j) + u[i][j];
            }
        }
    }

    public static void execM2() {
        // неявная схема
        prepareData();
        alpha[0] = 0.0;
        betta[0] = 0.0;

        for (int j = 0; j < m; j++) {
            for (int i = 1; i <= n; i++) {
                alpha[i] = -getBi(i-1) / (getAi(i-1) * alpha[i - 1] + getCi(i-1));
                betta[i] = (getFi(i-1, j) - getAi(i-1) * betta[i - 1]) / (getAi(i-1) * alpha[i - 1] + getCi(i-1));
            }

            u[n][j] = (getFi(n, j) - getAi(n) * betta[n]) / (alpha[n] * getAi(n) + getCi(n));

            for (int i = n - 1; i >= 0; i--) {
                u[i][j + 1] = u[i + 1][j] * alpha[i + 1] + betta[i + 1];
            }
        }
    }

    public static void execM3() {
        // схема Кранка-Никольсона
        prepareData();
        alpha[0] = 0.0;
        betta[0] = 0.0;
        for (int j = 0; j < m; j++) {
            for (int i = 1; i <= n; i++) {
                alpha[i] = -getBi_kn(i-1) / (getAi_kn(i-1) * alpha[i - 1] + getCi_kn(i-1));
                betta[i] = (getFi_kn(i-1, j) - getAi_kn(i-1) * betta[i - 1]) / (getAi_kn(i-1) * alpha[i - 1] + getCi_kn(i-1));
            }

            u[n][j] = (getFi_kn(n, j) - getAi_kn(n) * betta[n]) / (alpha[n] * getAi_kn(n) + getCi_kn(n));

            for (int i = n - 1; i >= 0; i--) {
                u[i][j + 1] = u[i + 1][j] * alpha[i + 1] + betta[i + 1];
            }
        }
    }

    private static void execE(String method) {
        double E = 0;
        double dif = 0;
        for (int j = 0; j <= m; j++) {
            for (int i = 0; i <= n; i++) {
                dif = Math.abs(getUxt(i, j) - u[i][j]);
                if (dif > E) {
                    E = dif;
                }
            }
        }
        System.out.println(method + ": h=" + h + " tau=" + tau + " E=" + E);
    }

    public static double getXi(int i) {
        return h * i;
    }

    public static double getTj(int j) {
        return tau * j;
    }

    public static double getMu(int j) {
        double t = getTj(j);
        return Math.log(t * t + 1);
    }

    public static double getTauGij(int i, int j) {
        double tj = getTj(j);
        return tau * (Math.sin(getXi(i)) + 2 * tj / (tj * tj + 1));
    }

    public static double getUxt(int i, int j) {
        return Math.sin(getXi(i) + getMu(j));
    }

    // ================================================
    public static Double getAi(int i) {
        return i == 0 || i == n ? 0 : A;
    }

    public static Double getBi(int i) {
        return i == 0 || i == n ? 0 : B;
    }

    public static Double getCi(int i) {
        return i == 0 || i == n ? 1 : C;
    }

    public static Double getFi(int i, int j) {
        if (i == 0 || i == n) {
            return getMu(j);
        }
        return -u[i][j] - getTauGij(i, j);
    }

    // ================================================
    public static Double getAi_kn(int i) {
        return i == 0 || i == n ? 0 : A_kn;
    }

    public static Double getBi_kn(int i) {
        return i == 0 || i == n ? 0 : B_kn;
    }

    public static Double getCi_kn(int i) {
        return i == 0 || i == n ? 1 : C_kn;
    }

    public static Double getFi_kn(int i, int j) {
        if (i == 0 || i == n) {
            return getMu(j);
        }
        return -(A_kn * u[i - 1][j] + c_kn * u[i][j] + B_kn * u[i + 1][j] + getTauGij(i, j));
    }

    // #############################################################
    public static void printU(Double[][] u) {
        for (int j = 0; j <= m; j++) {
            for (int i = 0; i <= n; i++) {
                System.out.print(u[i][j] + "  ");
            }
            System.out.println("");
        }
    }
}
