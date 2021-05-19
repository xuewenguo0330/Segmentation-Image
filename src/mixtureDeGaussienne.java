import javax.imageio.ImageIO;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Random;

public class mixtureDeGaussienne {

    /*calculer la distance entre un point d1 et le centre*/
    private static double dist(double[] d1, double[] centre) {
        double dist2 = 0;
        for(int i = 0; i < d1.length; i++) {
            dist2 += Math.pow((d1[i]-centre[i]), 2);
        }
        return Math.sqrt(dist2);
    }

    private static double distMin(double[] X, double[][] centres, int idx) {
        double dmin = dist(X,centres[0]);
        for(int c = 1; c < idx; c++) {
            double d = dist(X, centres[c]);
            if(d < dmin) { dmin = d;}
        }
        return dmin;
    }
    //I init
    public static double[][] initCentre(double[][] X , int nbrCentres) {
        double[][] centres = new double[nbrCentres][X[0].length];
        //initialise aléatoirement le 1er centre
        Random r1 = new Random();
        int idx0 = r1.nextInt(X.length-1);
        centres[0] = Arrays.copyOf(X[idx0], X[idx0].length);

        //pour les autres centres
        for(int c = 1; c < centres.length; c++) {
            double[] dist_min = new double[X.length];
            for(int i = 1; i < X.length; i++) {
                //la distance entre la donnee X[i] et la centre la plus proche
                dist_min[i] = distMin(X[i], centres, c);
            }
            double dmax = dist_min[0];
            int idx = 0;//l'indix du prochain centre
            for(int j = 1; j< dist_min.length; j++) {
                if(dmax < dist_min[j]) {
                    dmax = dist_min[j];
                    idx = j;
                }
            }
            centres[c] = Arrays.copyOf(X[idx], X[idx].length);
        }
        return centres;
    }

    /*II pour assignement*/
    private static double probaUnite(double[] X, double[] centre, double densite, double[] echelleCarre) {
        double res = densite;
        for(int i = 0; i< X.length; i++) {
            res *= (1./Math.sqrt(Math.PI*2*echelleCarre[i]))*Math.exp(Math.pow(X[i]-centre[i],2)/(-2*echelleCarre[i]));
        }
        //System.out.println((res));
        return res;
    }

    //k indice de centre (i.e de cluster)
    //Calculer la probabitité pour que x appartienne au cluster k
    private static double probaGauss(double[] X, double[][] centres, double[] densites, double[][] echellesCarre, int k) {
        double proba = probaUnite(X, centres[k], densites[k], echellesCarre[k]);
        //System.out.println("hello" + densites[k]);
        double somme = proba;

        for(int ctr = 0; ctr < centres.length; ctr++ ) {
            if(ctr!=k) {
                somme += probaUnite(X, centres[ctr], densites[ctr], echellesCarre[ctr]);
            }
        }
        //System.out.println("hello " + proba/somme);
        return proba/somme;
    }

    /*Assigner
     * X        une liste de points dans un espace à D dimensions
     * centres  chaque centre (parmi les k centres) par une position dans l'espace à D dimensions
     */
    public static double[][] Assigner(double[][] X, double[][] centres, double[] densites, double[][] echellesCarre) {
        double[][] ass = new double[X.length][centres.length];
        //pour chaque pts dans X
        for(int i = 0; i< X.length; i++) {
            for(int c =0; c< centres.length; c++) {
                ass[i][c] = probaGauss(X[i], centres, densites, echellesCarre, c);
            }
        }
        //test
        for(int i = 0; i< X.length; i++) {
            double s = 0;
            for(int c =0; c< centres.length; c++) {
                s += ass[i][c];
            }
            //System.out.println("s= " + s +" ass[0] ="+ Arrays.toString(ass[X.length-5]));
            //somme de probabilite varie entre 0.99999999, 1 et 1.0000000000000002
            if(Math.abs(s-1)>0.000000000001){
                System.out.println("erreur somme proba != 1. indice =" + i + " s = " + s);
            }
        }
        return ass;
    }

    /*III Deplct et maj*/
    private static double[] Rk(double[][] ass) {
        double[] Rk = new double[ass[0].length];
        for (int c = 0; c < ass[0].length; c++) {
            Rk[c] = 0;
            for (int i = 0; i < ass.length; i++) {
                Rk[c] += ass[i][c];
            }
        }
        return Rk;
    }

    private static double majCentre(double[][] X, double[][] centres, double[][] ass, double[] RK) {
        double distDep = 0;
        for(int c =0; c < centres.length; c++) {//pour chaque cluster(centre)
            double[] ancienCtr = Arrays.copyOf(centres[c], centres[c].length);

            for(int j =0; j < centres[0].length; j++) {//pour chaque coord du centre
                double coordCentre_cj = 0;
                for(int i =0; i < X.length; i++){//on rajoute tous les rc*Xi,j
                    coordCentre_cj += X[i][j]*ass[i][c];
                }
                centres[c][j] = coordCentre_cj/RK[c];//on met a jour les centres
            }
            distDep += dist(ancienCtr,centres[c]);
        }
        return distDep;
    }

    private static void majEcheCarre(double[][] X, double[][] centres, double[][] ass, double[][] echelelsCarre, double[] RK) {
        for(int c =0; c < centres.length; c++) {//pour chaque cluster(centre)
            for(int j =0; j < centres[0].length; j++) {//pour chaque coord du centre

                double coordEchelle_cj = 0;
                for(int i =0; i < X.length; i++){//on rajoute tous les (ri,c)*(Xi,j-Mc,j)
                    coordEchelle_cj += Math.pow((X[i][j]-centres[c][j]),2)*ass[i][c];
                }
                echelelsCarre[c][j] = coordEchelle_cj/RK[c];//on met a jour les echelles
            }
        }
    }

    private static void majDensite(double[][] X, double[] RK, double[] densites){
        for(int i = 0; i< RK.length; i++) {
            densites[i] = densites[i]/X[0].length;
        }
    }

    public static double Deplct(double[][] X, double[][] centres, double[][] ass,double[] densite,double[][] echellesCarre) {
        double[] RK = Rk(ass);
        double distDep = majCentre(X,centres,ass,RK);
        //System.out.println("distDep = " + distDep);
        majEcheCarre(X,centres,ass, echellesCarre,RK);
        majDensite(X,RK,densite);
        return distDep;
    }

    public static int maxProbaIdx(double[] data_ass){
        double proba = data_ass[0];
        int res = 0;
        for(int i = 1; i< data_ass.length; i++){
            if(proba < data_ass[i]){
                proba = data_ass[i];
                res = i;
            }
        }
        return res;
    }
    public static double[][] loadImgData(BufferedImage bui) throws IOException {
        int width = bui.getWidth();
        int height = bui.getHeight();
        System.out.println("Hauteur=" + width);
        System.out.println("Largeur=" + height);

        double[][] rgbData = new double[width*height][3];
        for(int i = 0; i<width; i++){
            for(int j = 0; j<height; j++){
                Color c = new Color(bui.getRGB(i, j));
                // Calcul des trois composant de couleurs inverse et normalisé à 1
                rgbData[i+j*width][0] = (double) (255-c.getRed())/255.0;
                rgbData[i+j*width][1] = (double) (255-c.getGreen())/255.0;
                rgbData[i+j*width][2]  = (double) (255-c.getBlue())/255.0;
            }
        }
        System.out.println("RGB normalisé= "+ Arrays.toString(rgbData[0]));
        return  rgbData;
    }

    public static void main(String[] args) throws IOException {
        // charger les donnees
        String path = "./";
        String imageMMS = path + "mms.png";
        BufferedImage bui = ImageIO.read(new File(imageMMS));

        int width = bui.getWidth();
        int height = bui.getHeight();
        double[][] rgbData =loadImgData(bui); // width*hight points en 3D


        //initialiser les centres, les densites et les echelles
        int D = 3; // 3 dimensions = R G B
        int k = 8; // 8 centres
        //double[][] centres = initCentre(rgbData,k);

        //initialiser les centres pas aleatoirement pour l'image inverse
        double[][] centres = new double[k][D];
        centres[0][0] = 1;   centres[0][1] = 1;  centres[0][2] = 1; //blanc inverse de noir
        centres[1][0] = 0;   centres[1][1] = 1; centres[1][2] = 1; //cyan inverse de rouge
        centres[2][0] = 1;   centres[2][1] = 0; centres[2][2] = 1; //magenta inverse de vert
        centres[3][0] = 1;   centres[3][1] = 1; centres[3][2] = 0; //jaune inverse de bleu
        centres[4][0] = 0;  centres[4][1] = 0; centres[4][2] = 1;  //bleu inverse de jaune
        centres[5][0] = 0;  centres[5][1] = 1;centres[5][2] = 0; //vert inverse pour un autre rouge dans mms
        centres[6] = rgbData[0]; //la couleur du fond
        centres[7][0] = 0;  centres[7][1] = 0;centres[7][2] = 0; //noir inverse du blanc pour le reflet

        for(int i = 0; i< centres.length; i++){
            System.out.println("centresINit["+ i + "] = "+Arrays.toString((centres[i])));
        }

        double[] densites = new double[k];
        double[][] echellesCarre = new double[k][D];
        for(int i = 0; i < k; i++){ //initialiser tous les densite à 1./k
            densites[i] = 1./k;
            for(int j=0; j< D; j++){//initialiser tous les echelle à 1
                echellesCarre[i][j] = 1;
            }
        }

        //la mixture de Gaussienne
        double eps=0.001;
        double maj = 10;
        double[][] ass = new double[rgbData.length][centres.length];
        while(maj>eps) {
            ass = Assigner(rgbData,centres,densites,echellesCarre);
            maj = Deplct(rgbData,centres,ass,densites,echellesCarre);
        }

        // verification
        for(int i=0; i<centres.length; i++) {
            System.out.println("Pos centre "+i+": "+centres[i][0]+" "+centres[i][1]+ " "+centres[i][2]);
        }

        //test Image
        Color[] centersColor = new Color[k];
        Color[] tabColor = new Color[rgbData.length];

        for (int i = 0; i < k ; i++) {
            centersColor[i] = new Color((float) centres[i][0], (float) centres[i][1], (float) centres[i][2]);
            System.out.println("Pos centreColor "+i+": "+centersColor[i].toString());
        }
        for (int i = 0; i < rgbData.length ; i++) {
            tabColor[i] = centersColor[maxProbaIdx(ass[i])];
        }

        BufferedImage bui_out = new BufferedImage(bui.getWidth(),bui.getHeight(),BufferedImage.TYPE_3BYTE_BGR);
        BufferedImage bui_out0 = new BufferedImage(bui.getWidth(),bui.getHeight(),BufferedImage.TYPE_3BYTE_BGR);
        for(int i=0 ; i<height ; i++) {
            for(int j=0 ; j<width ; j++) {
                bui_out.setRGB(j, i, tabColor[i * width + j].getRGB());
                bui_out0.setRGB(j, i, 255-tabColor[i * width + j].getRGB());
            }
        }

        ImageIO.write(bui_out, "PNG", new File(path+"InverseParti1.png"));
        ImageIO.write(bui_out0, "PNG", new File(path+"Parti1.png"));
    }
}
