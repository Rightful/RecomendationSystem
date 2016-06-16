import java.util.*;

public class main {

    public static void main(String[] args) {
        Locale.setDefault(Locale.US);
        // Create user list, movie list, and list of ratings
        UserList userList = new UserList();
        userList.readFile("data/users.csv");
        MovieList movieList = new MovieList();
        movieList.readFile("data/movies.csv");
        RatingList ratingList = new RatingList();
        ratingList.readFile("data/ratings.csv", userList, movieList);
        // Read list of ratings we need to predict
        RatingList predRatings = new RatingList();
        predRatings.readFile("data/predictions.csv", userList, movieList);

        user_similar(userList, movieList, ratingList, predRatings);

        double item[] = item_based(userList, movieList, ratingList, predRatings);
       // double user[] = user_based(userList, movieList, ratingList, predRatings);

        for (int i = 0; i < predRatings.size(); i++) {
            double itemm = item[i];
            //double userr = user[i];
//            double item_user = 0.1 * userr + 0.9 * itemm;
            double item_user = itemm;
            predRatings.get(i).setRating(item_user);
        }


        predRatings.writeResultsFile("submission.csv");
    }

    private static double mean_all = 2.4513148902975604;

    public static double[][] movie_user(UserList userList,
                                        MovieList movieList, RatingList ratingList, RatingList predRatings) {

        double[][] movie_user = new double[movieList.size() + 1][userList.size() + 1];


        for (int z = 0; z < ratingList.size(); z++) {
            movie_user[ratingList.get(z).getMovie().getIndex()][ratingList.get(z).getUser().getIndex()] = ratingList.get(z).getRating();
        }
        return movie_user;
    }

    public static double[][] rating_deviation_m_u(UserList userList,
                                                  MovieList movieList, RatingList ratingList, RatingList predRatings) {

        double[][] m_u = movie_user(userList, movieList, ratingList, predRatings);
        double[][] deviation_m_u = new double[movieList.size() + 1][userList.size() + 1];
        double[] user = new double[userList.size() + 1];
        double[] movie = new double[movieList.size() + 1];

        double sum_movie = 0;
        double sizem = 0;
        for (int i = 1; i < movieList.size() + 1; i++) {
            sum_movie = 0;
            sizem = 0;
            for (int j = 1; j < userList.size() + 1; j++) {
                sum_movie += m_u[i][j];
                if (m_u[i][j] != 0.0) {
                    sizem++;
                }
            }
            movie[i] = (sum_movie / sizem) - mean_all;
        }

        double sum_user = 0;
        double sizeu = 0;
        for (int i = 1; i < userList.size() + 1; i++) {
            sum_user = 0;
            sizeu = 0;
            for (int j = 1; j < movieList.size() + 1; j++) {
                sum_user += m_u[j][i];
                if (m_u[j][i] != 0.0) {
                    sizeu++;
                }
            }
            user[i] = (sum_user / sizeu) - mean_all;
        }


        for (int j = 1; j < movieList.size() + 1; j++) {
            deviation_m_u[j][0] = movie[j];
        }

        for (int j = 1; j < userList.size() + 1; j++) {
            deviation_m_u[0][j] = user[j];
        }
        return deviation_m_u;
    }

    public static double[][] meaned_movie_user(UserList userList,
                                               MovieList movieList, RatingList ratingList, RatingList predRatings) {

        double[][] meaned_movie_user = movie_user(userList, movieList, ratingList, predRatings);
        double[] mean_array = new double[movieList.size() + 1];

        fillMeanArray(userList, movieList, meaned_movie_user, mean_array);
        fillMovieUser(userList, movieList, meaned_movie_user, mean_array);

        return meaned_movie_user;
    }

    private static void fillMovieUser(UserList userList, MovieList movieList, double[][] meaned_movie_user, double[] mean_array) {
        double rating;
        for (int m = 1; m < movieList.size() + 1; m++) {
            for (int u = 1; u < userList.size() + 1; u++) {
                rating = meaned_movie_user[m][u];
                if (meaned_movie_user[m][u] == 0.0) {
                    meaned_movie_user[m][u] = 0.0;
                } else {
                    meaned_movie_user[m][u] = (rating - mean_array[m]);
                }
            }
        }
    }

    private static void fillMeanArray(UserList userList, MovieList movieList, double[][] meaned_movie_user, double[] mean_array) {
        double sum;
        double size;
        double mean;
        for (int m = 1; m < movieList.size() + 1; m++) {
            sum = 0;
            size = 0;
            for (int u = 1; u < userList.size(); u++) {
                sum += meaned_movie_user[m][u];
                if (meaned_movie_user[m][u] > 0) {
                    size++;
                }
            }
            mean = sum / size;
            mean_array[m] = mean;
        }
    }

    public static double[][] user_similar(UserList userList,
                                          MovieList movieList, RatingList ratingList, RatingList predRatings) {

        double user_similar[][] = meaned_movie_user(userList, movieList, ratingList, predRatings);
        double sim[][] = new double[userList.size() + 1][userList.size() + 1];

        double dotproduct = 0;
        double inorm = 0;
        double jnorm = 0;
        double simij = 0;
        for (int s = 1; s < userList.size() + 1; s++) {
            System.out.print("\r" + s);
            for (int m = 2; m < userList.size() + 1; m++) {
                dotproduct = 0;
                inorm = 0;
                jnorm = 0;
                simij = computeSimilarity(movieList, user_similar, dotproduct, inorm, jnorm, s, m);
                if (simij >= 0) {
                    sim[s][m] = Double.isNaN(simij) ? -1 : simij;
                } else {
                    sim[s][m] = Double.isNaN(simij) ? -1 : -1;
                }
            }
        }

        // similarity van zelfde user row is altijd 1
        // maar deze wil je niet meenemen dus op -1
        for (int m = 1; m < userList.size() + 1; m++) {
            sim[m][m] = -1;
        }
        return sim;
    }

    private static double computeSimilarity(MovieList movieList, double[][] user_similar, double dotproduct, double inorm, double jnorm, int s, int m) {
        double i;
        double j;
        double simij;
        for (int u = 1; u < movieList.size() + 1; u++) {
            i = user_similar[u][s];
            j = user_similar[u][m];
            dotproduct += i * j;
            inorm += Math.pow(i, 2);
            jnorm += Math.pow(j, 2);
        }
        simij = dotproduct / (Math.sqrt(inorm) * Math.sqrt(jnorm));
        return simij;
    }

    public static double[][] movie_similar(UserList userList,
                                           MovieList movieList, RatingList ratingList, RatingList predRatings) {

        double movie_similar[][] = meaned_movie_user(userList, movieList, ratingList, predRatings);
        double sim[][] = new double[movieList.size() + 1][movieList.size() + 1];

        double i = 0;
        double j = 0;
        double dotproduct = 0;
        double inorm = 0;
        double jnorm = 0;
        double simij = 0;

        for (int s = 1; s < movieList.size() + 1; s++) {
            System.out.print("\r" + s);
            for (int m = 2; m < movieList.size() + 1; m++) {
                dotproduct = 0;
                inorm = 0;
                jnorm = 0;
                simij = computeSimilarity2(userList, dotproduct, inorm, jnorm, movie_similar[s], movie_similar[m]);
                if (simij >= 0) {
                    sim[s][m] = Double.isNaN(simij) ? -1 : simij;
                } else {
                    sim[s][m] = Double.isNaN(simij) ? -1 : -1;
                }
            }
        }

        // similarity van zelfde movie row is altijd 1
        // maar deze wil je niet meenemen dus op -1
        for (int m = 1; m < movieList.size() + 1; m++) {
            sim[m][m] = -1;
        }
        return sim;
    }

    private static double computeSimilarity2(UserList userList, double dotproduct, double inorm, double jnorm, double[] doubles, double[] doubles1) {
        double i;
        double j;
        double simij;
        for (int u = 1; u < userList.size() + 1; u++) {
            i = doubles[u];
            j = doubles1[u];
            dotproduct += i * j;
            inorm += Math.pow(i, 2);
            jnorm += Math.pow(j, 2);
        }
        simij = dotproduct / (Math.sqrt(inorm) * Math.sqrt(jnorm));
        return simij;
    }

    public static double[] item_based(UserList userList,
                                      MovieList movieList, RatingList ratingList, RatingList predRatings) {

        double matrix[][] = movie_user(userList, movieList, ratingList, predRatings);
        double sim[][] = movie_similar(userList, movieList, ratingList, predRatings);
        double rating_array[] = new double[predRatings.size()];
        double dev_m_u[][] = rating_deviation_m_u(userList, movieList, ratingList, predRatings);

        HashMap<Integer, TreeMap<Double, Integer>> movie_i_sim_movie_j = new HashMap<>();

        double similarity_value = 0;
        for (int j = 0; j < predRatings.size(); j++) {
            TreeMap<Double, Integer> sim_movie = new TreeMap<>(Collections.reverseOrder());
            for (int i = 1; i < movieList.size() + 1; i++) {
                similarity_value = sim[predRatings.get(j).getMovie().getIndex()][i];
                sim_movie.put(similarity_value, i);
            }
            movie_i_sim_movie_j.put(j, sim_movie);

            int reached = 0;
            int k = 18;
            double rating = 0;
            double ratingsimcount = 0;
            double simcount = 0;
            double rxj = 0;
            double sij = 0;
            double bxi = mean_all + dev_m_u[0][predRatings.get(j).getUser().getIndex()] + dev_m_u[predRatings.get(j).getMovie().getIndex()][0];
            for (Map.Entry<Double, Integer> entry : sim_movie.entrySet()) {

                if (reached == k) {
                    break;
                }

                if (matrix[entry.getValue()][predRatings.get(j).getUser().getIndex()] == 0.0) {
                    continue;
                }
                double bxj = (mean_all + dev_m_u[0][predRatings.get(j).getUser().getIndex()] + dev_m_u[entry.getValue()][0]);
                rxj = matrix[entry.getValue()][predRatings.get(j).getUser().getIndex()];
                sij = entry.getKey();

                ratingsimcount += sij * (rxj - bxj);
                simcount += sij;
                reached++;
            }
            rating = bxi + ratingsimcount / simcount;

            if (rating > 5)
                rating = 5;
            if (rating < 1)
                rating = 1;

            rating_array[j] = Double.isNaN(rating) ? mean_all : rating;
            movie_i_sim_movie_j.remove(j); // free unused memory
        }
        return rating_array;
    }

    public static double[] user_based(UserList userList,
                                      MovieList movieList, RatingList ratingList, RatingList predRatings) {

        double matrix[][] = movie_user(userList, movieList, ratingList, predRatings);
        double sim[][] = user_similar(userList, movieList, ratingList, predRatings);
        double rating_array[] = new double[predRatings.size()];
        double dev_m_u[][] = rating_deviation_m_u(userList, movieList, ratingList, predRatings);

        HashMap<Integer, TreeMap<Double, Integer>> user_i_user_j = new HashMap<>();

        double similarity_value = 0;
        for (int j = 0; j < predRatings.size(); j++) {
            TreeMap<Double, Integer> sim_user = new TreeMap<>(Collections.reverseOrder());
            for (int i = 1; i < userList.size() + 1; i++) {
                similarity_value = sim[predRatings.get(j).getUser().getIndex()][i];
                sim_user.put(similarity_value, i);
            }
            user_i_user_j.put(j, sim_user);

            int reached = 0;
            int k = 18;
            double rating = 0;
            double ratingsimcount = 0;
            double simcount = 0;
            double rxj = 0;
            double sij = 0;
            double bxi = mean_all + dev_m_u[0][predRatings.get(j).getUser().getIndex()] + dev_m_u[predRatings.get(j).getMovie().getIndex()][0];
            for (Map.Entry<Double, Integer> entry : sim_user.entrySet()) {

                if (reached == k) {
                    break;
                }

                if (matrix[predRatings.get(j).getMovie().getIndex()][entry.getValue()] == 0.0) {
                    continue;
                }
                double bxj = (mean_all + dev_m_u[predRatings.get(j).getMovie().getIndex()][0] + dev_m_u[0][entry.getValue()]);
                rxj = matrix[predRatings.get(j).getMovie().getIndex()][entry.getValue()];
                sij = entry.getKey();

                ratingsimcount += sij * (rxj - bxj);
                simcount += sij;
                reached++;
            }
            rating = bxi + ratingsimcount / simcount;

            if (rating > 5)
                rating = 5;
            if (rating < 1)
                rating = 1;

            rating_array[j] = Double.isNaN(rating) ? mean_all : rating;
            user_i_user_j.remove(j); // free unused memory
        }
        return rating_array;
    }


}