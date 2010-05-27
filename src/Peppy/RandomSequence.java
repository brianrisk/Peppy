package Peppy;

import java.util.Random;

public class RandomSequence {
	
	public static void main(String args[]) {
		Random random = new Random();
		int choice;
		char out = 'A';
		for (int i = 0; i < 50; i++) {
			choice = random.nextInt(4);
			if (choice == 0) out = 'A';
			if (choice == 1) out = 'T';
			if (choice == 2) out = 'G';
			if (choice == 3) out = 'C';
			System.out.print(out);
		}
	}

}
