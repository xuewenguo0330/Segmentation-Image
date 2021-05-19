import javax.imageio.ImageIO;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;

public class SixMeans {


    private static double distance(double[] p1,double[] p2) {
        float somme=0;
        for(int i=0;i<p1.length;i++) {
            somme+=Math.pow(p1[i]-p2[i],2);
        }
        return Math.sqrt(somme);
    }
    public static int[]Assigner(double[][]X,double[][]centres){
        //X une liste de points dans un espace à D dimensions
        //centres une liste de centres parmi les k centres par une positon dans l'espace à D dimensions

        int[] ass = new int[X.length];
        //pour chaque pts dans X
        for(int i = 0; i< X.length; i++) {
            ass[i] = 0; // pour commencer on lui assigne le centre 0
            double distMin = distance(X[i], centres[0]); // on calcule la distance entre cette pts et c0
            //on cherche le centre le plus proche de cette points  et on lui assigne à ce centre
            for(int c = 1; c < centres.length; c++) {
                if(distance(X[i], centres[c])< distMin) {
                    ass[i] = c;
                }
            }
        }
        return ass;
    }
    private static double[]NouvelleCentre(double [][]X,int[]ass,int l) {//ass, 每个点属于哪个centre
        double[] nouvelle = new double[X[0].length];
        float sommex = 0.f;
        float sommey = 0.f;
        float fois = 0.f;
        for (int i = 0; i < X.length; i++) {
            if (ass[i] == l) {
                sommex += X[i][0];
                sommey += X[i][1];
                fois++;
            }
        }
        nouvelle[0] = sommex / fois;
        nouvelle[1] = sommey / fois;
        return nouvelle;
    }


    public static double Deplct(double[][]X,double centres[][],int[]ass) {
        double somme=0;//somme de distances
        double d;
        for(int i=0;i<centres.length;i++) {
            double[] NC=NouvelleCentre(X,ass,i);
            d=distance(centres[i],NC);
            somme+=d;
            centres[i]=NC;
        }

        return somme/centres.length;
    }




    public static void main(String[] args)  throws IOException {
        String path = "./";
        String imageMMS = path + "mms.png";

        // Lecture de l'image ici
        BufferedImage bui = ImageIO.read(new File(imageMMS));//read image
        int width = bui.getWidth();
        int height = bui.getHeight();

        int[] im_pixels = bui.getRGB(0, 0, width, height, null, 0, width);
        int M=im_pixels.length; //nombre de pixels, et donc taille du tableau image

        int D=3; // trois dimensions
        int k=8; // huit centres(6 couleurs+fond+reflet )
        double[][] X = new double[M+1][D]; // M points en D dimensions
        double[][] centres = new double[k][D];

        //creation des tableau de pixel
        System.out.println("M "+ M);
        for(int i=0;i<M;i++){
            Color c=new Color(im_pixels[i]);
            X[M][0]=(double) c.getRed()/255;
            X[M][1]=(double) c.getGreen()/255;
            X[M][2]=(double) c.getBlue()/255;
        }

        //creation des centre
        Color c = new Color(im_pixels[0]);         //couleur du fond
        centres[0][0]=0.f;  centres[0][1]=0.f;  centres[0][2]=1.f;    //bleu
        centres[1][0]=1.f;  centres[1][1]=1.f;  centres[1][2]=0.f;    //jaune
        centres[2][0]=0.f;  centres[2][1]=1.f;  centres[2][2]=0.f;    //vert
        centres[3][0]=1.f;  centres[3][1]=0.5;  centres[3][2]=0.f;  //orange
        centres[4][0]=0.f;  centres[4][1]=0.f;  centres[4][2]=0.f;    //noir
        centres[5][0]=1.f;  centres[5][1]=0.f;  centres[5][2]=0.f;    //rouge
        centres[6][0]=(double) c.getRed()/255.0;  centres[6][1]=(double) c.getGreen()/255.0;  centres[6][2]=(double) c.getBlue()/255.0;    //fond
        centres[7][0]=1.f;  centres[7][1]=1.f;  centres[7][2]=1.f;    //reflets

        //segement par couleur
        double eps=0.001;
        double maj = 10;
        int[] ass = new int[X.length];
        while(maj>eps) {
            ass = Assigner(X,centres);
            maj = Deplct(X,centres,ass);
            eps=10;
        }

        // verification
        for(int i=0; i<X.length; i++) {
            System.out.println("Pt "+i+" assigné à "+ass[i]);
        }
        for(int i=0; i<centres.length; i++) {
            System.out.println("Pos centre "+i+": "+centres[i][0]+" "+centres[i][1]);
        }
    }

}
