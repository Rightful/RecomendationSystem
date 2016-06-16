import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

public class RatingList extends ArrayList<Rating> {

	// Reads in a file with rating data
	public void readFile(String filename, UserList userList, MovieList movieList) {
		BufferedReader br = null;
		String line;
		try {
			br = new BufferedReader(new FileReader(filename));
			while ((line = br.readLine()) != null) {
				String[] ratingData = line.split(";");
				if (ratingData.length == 3) {
					add(new Rating(
							userList.get(Integer.parseInt(ratingData[0]) - 1),
							movieList.get(Integer.parseInt(ratingData[1]) - 1),
							Double.parseDouble(ratingData[2])));
				} else {
					add(new Rating(
							userList.get(Integer.parseInt(ratingData[0]) - 1),
							movieList.get(Integer.parseInt(ratingData[1]) - 1),
							0.0));
				}
			}
		} catch (Exception e) {
			e.printStackTrace();
		} finally {
			if (br != null) {
				try {
					br.close();
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		}
	}

	// Writes a result file
	public void writeResultsFile(String filename) {
		PrintWriter pw;
		try {
			pw = new PrintWriter(filename);
			pw.println("Id,Rating");
			for (int i = 0; i < size(); i++) {
				pw.println((i + 1) + "," + get(i).getRating());
			}
			pw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}
}
