import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;

public class UserList extends ArrayList<User> {

	// Reads in a file with user data
	public void readFile(String filename) {
		BufferedReader br = null;
		String line;
		try {
			br = new BufferedReader(new FileReader(filename));
			while ((line = br.readLine()) != null) {
				String[] userData = line.split(";");
				add(Integer.parseInt(userData[0]) - 1,
						new User(Integer.parseInt(userData[0]), userData[1]
								.equals("M"), Integer.parseInt(userData[2]),
								Integer.parseInt(userData[3])));
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
}
