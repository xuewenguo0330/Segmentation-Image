
import java.util.ArrayList;
import javax.imageio.ImageIO;
import java.awt.image.BufferedImage;
import java.awt.Color;
import java.io.File;
import java.io.IOException;

public class MixGauss {



	public static double fonction_gaussienne(double x, double m, double sigma) {
		return Math.exp(-Math.pow(x-m,2)/(2*Math.pow(sigma, 2)))/Math.sqrt(2*Math.PI*Math.pow(sigma, 2));
	}


	public static double[][] maj_m(double[][] m, double[][] X, double[][] probas) {
		int K = probas[0].length;
		int D = X[0].length;
		int M = X.length;
		for (int k =0; k<K; k++) {
			double R=0;
			for (int d = 0; d<M; d++) {
				R+=probas[d][k];
			}
			for (int i=0; i<D; i++)  {
				double S = 0;
				for (int d=0; d<M; d++) {
					S+=probas[d][k]*X[d][i];
				}
				m[k][i]=S/R;
			}
		}
		return m;
	}


	public static double[][] maj_sigma(double[][] sigma, double[][] X, double[][] probas, double[][] m) {
		int K = probas[0].length;
		int D = X[0].length;
		int M = X.length;
		for (int k =0; k<K; k++) {
			double R=0;
			for (int d = 0; d<M; d++) {
				R+=probas[d][k];
			}
			for (int i=0; i<D; i++)  {
				double S = 0;
				for (int d=0; d<M; d++) {
					S+=probas[d][k]*Math.pow(X[d][i]-m[k][i],2);
				}
				sigma[k][i]=Math.sqrt(S/R);
			}
		}
		return sigma;
	}


	public static double[] maj_rho(double[] rho, double[][] X, double[][] probas) {
		int K = probas[0].length;
		int M = X.length;
		for (int k =0; k<K; k++) {
			double R=0;
			for (int d = 0; d<M; d++) {
				R+=probas[d][k];
			}
			rho[k]=R/M;
		}
		return rho;
	}


	public static int Indice_max(double[] tab) {
		int imax=0;
		for (int i=1; i<tab.length; i++) {
			if (tab[i]>tab[imax]) {
				imax=i;
			}
		}
		return imax;
	}


	public static double Distance(double[] X, double[] Y) {
		double D = 0;
		for (int i=0; i<X.length;i++) {
			D+=Math.pow(X[i]-Y[i],2);
		}
		return Math.sqrt(D);
	}


	public static void maj_probas(double[][] X, double[][] probas, double[][] m, double[][] sigma, double[] rho) {
		int M= X.length;             //nombre de donn?es
		int K = probas[0].length;    //nombre de centres ("probas" correspond au "Assigner()" pr?c?dent)
		int D = X[0].length;         //dimension des donn?es
		for (int d=0; d<M; d++) {
			double S = 0;
			for (int l=0; l<K; l++) {
				double G = 1;
				for (int i=0; i<D; i++) {
					G*= fonction_gaussienne(X[d][i], m[l][i], sigma[l][i]);
				}
				G*=rho[l];
				S+=G;
			}
			for (int k=0; k<K; k++) {
				double P = 1;
				for (int i=0; i<D; i++) {
					P*= fonction_gaussienne(X[d][i], m[k][i], sigma[k][i]);
				}
				P*=rho[k];
				probas[d][k]=P/S;
			}
		}
	}


	public static double Deplct(double[][] X, double[][] probas, double[][] m, double[][] sigma, double[] rho) {
		int K = probas[0].length;
		int D = X[0].length;
		double[][] m0 = new double[K][D];
		for (int k=0; k<K; k++) {
			for (int i=0; i<D; i++) {
				m0[k][i]=m[k][i];
			}
		}
		maj_m(m,X,probas);
		maj_sigma(sigma,X,probas,m);
		maj_rho(rho,X,probas);
		double dist = 0;
		for (int k=0; k<K; k++) {
			dist+=Distance(m0[k],m[k]);
		}
		return dist;
	}


	static BufferedImage bi;

	public static void main(String[] args) throws IOException {
		int D=3; // trois dimensions
		int K=3; // huit centres (6 couleurs + fond + reflets blanc)


		String path = "./";
		String imageMMS = path + "mms.png";

		// Lecture de l'image ici
		BufferedImage bui = ImageIO.read(new File(imageMMS));

		int width = bui.getWidth();
		int height = bui.getHeight();

		int[] im_pixels = bui.getRGB(0, 0, width, height, null, 0, width);
		int M=im_pixels.length; //nombre de pixels, et donc taille du tableau image

		/**cr?ation du tableau de pixels**/
		double[][] image = new double[M][D];   // tableau des pixels en couleur 3D
		for(int i=0 ; i<M ; i++) {
			Color c = new Color(im_pixels[i]);
			image[i][0] = (double) c.getRed()/255.0;
			image[i][1] = (double) c.getGreen()/255.0;
			image[i][2] = (double) c.getBlue()/255.0;
		}


		/**cr?ation des centres**/
		double[][] m = new double[K][D];
		Color c = new Color(im_pixels[0]);         //couleur du fond

		m[0][0]=0;  m[0][1]=0;  m[0][2]=1;    //bleu

		m[1][0]=1;  m[1][1]=1;  m[1][2]=0;    //jaune
		m[2][0]=0;  m[2][1]=1;  m[2][2]=0;    //vert
		m[3][0]=1;  m[3][1]=0.5;  m[3][2]=0;  //orange
		m[4][0]=0;  m[4][1]=0;  m[4][2]=0;    //noir
		m[5][0]=1;  m[5][1]=0; m[5][2]=0;    //rouge
		m[6][0]=(double) c.getRed()/255.0;  m[6][1]=(double) c.getGreen()/255.0;  m[6][2]=(double) c.getBlue()/255.0;    //fond
		m[7][0]=1;  m[7][1]=1;  m[7][2]=1;    //reflets


		/**initialisation des param?tres des gaussiennes**/
		double[][] sigma = new double[K][D];
		double[] rho = new double[K];
		double[][] probas = new double[M][K];

		for (int d=0; d<M; d++) {
			for (int k=0; k<K; k++) {
				probas[d][k]=1/K;
			}
		}


		for (int k=0; k<K; k++) {
			rho[k]=(double) 1/K;
			System.out.println(rho[k]);
		}

		for (int k=0; k<K; k++) {
			for (int i=0; i<D; i++) {
				if(k<6) {      //couleur
					sigma[k][i]=Math.sqrt(0.3);
				}
				else {         //fond ou reflets
					sigma[k][i]=Math.sqrt(0.5);
				}
			}
		}


		/**segmentation par couleur**/
		double eps=0.001;
		double maj = 10;

		while(maj>eps) {
			maj_probas(image,probas,m,sigma,rho);
			maj = Deplct(image,probas,m,sigma,rho);
		}

		/**score**/
		double[] score = new double[M];
		for (int d=0; d<M; d++) {
			double S=0;
			for (int k=0; k<K; k++) {
				double P=1;
				for (int i=0; i<D; i++) {
					P*=fonction_gaussienne(image[d][i],m[k][i],sigma[k][i]);
				}
				S+=rho[k]*P;
			}
			score[d]=Math.log(S);
		}

		double S=0;
		for (int d=0; d<M; d++) {
			S+=score[d];
		}
		double L = S/M;


		/**verification**/
		System.out.println("score :"+L);         //score total
		for (int k=0; k<K; k++) {
			//System.out.println(rho[k]);
		}

		Color[] tabColor= new Color[M];
		for(int d=0 ; d<M ; d++) {
			double[] BonCentre = m[Indice_max(probas[d])];
			tabColor[d]=new Color((int)(BonCentre[0]*255),(int)(BonCentre[1]*255),(int)(BonCentre[2]*255));
		}

		BufferedImage bui_out = new BufferedImage(bui.getWidth(),bui.getHeight(),BufferedImage.TYPE_3BYTE_BGR);
		for(int i=0 ; i<height ; i++) {
			for(int j=0 ; j<width ; j++) {
				bui_out.setRGB(j,i,tabColor[i*width+j].getRGB());
			}
		}
		ImageIO.write(bui_out, "PNG", new File(path+"test.png"));         //affichage de l'image constitu?e des pixels de la couleur la plus probable d'apr?s l'algorithme

	}
}
