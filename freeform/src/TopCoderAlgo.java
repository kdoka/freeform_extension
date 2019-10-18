import java.util.Arrays;


public class TopCoderAlgo {
	
	
	
	public static int INF = 100000000 ;   //just infinity
	int n;
	double[] lx; 
	double[] ly;        		//labels of X and Y parts
	int[] xy;               //xy[x] - vertex that is matched with x,
	int[] yx;               //yx[y] - vertex that is matched with y
	boolean[] S;
	boolean[] T;		//sets S and T in algorithm
	double[] slack;            //as in the algorithm description
	int[] slackx;           //slackx[y] such a vertex, that
	                         // l(slackx[y]) + l(y) - w(slackx[y],y) = slack[y]
	int[] prev;             //array for memorizing alternating paths
	
	int ret;                      //weight of the optimal matching
    int max_match;                    //number of vertices in current matching
		
	void hungarian(double[][] cost, int[] assignment){
		n = cost.length;
		//max_match;        //n workers and n jobs
		lx = new double[n]; 
		ly = new double[n];        		//labels of X and Y parts
		xy = new int[n];               //xy[x] - vertex that is matched with x,
		yx = new int[n];               //yx[y] - vertex that is matched with y
		S = new boolean[n];
		T = new boolean[n];		//sets S and T in algorithm
		slack = new double[n];            //as in the algorithm description
		slackx = new int[n];           //slackx[y] such a vertex, that
		                         // l(slackx[y]) + l(y) - w(slackx[y],y) = slack[y]
		prev=new int[n];             //array for memorizing alternating paths
		
		ret = 0;                      //weight of the optimal matching
	    max_match = 0;                    //number of vertices in current matching
	    
	    Arrays.fill(xy, -1);
	    Arrays.fill(yx, -1);
	    
	   //init labels 
	    Arrays.fill(lx, 0.0);
	    Arrays.fill(ly, 0.0);
	    for (int x = 0; x < n; x++)
	        for (int y = 0; y < n; y++)
	            lx[x] = Math.max(lx[x], -cost[x][y]);                   //step 0
	   
	    augment(cost); 
	    //steps 1-3
	    for (int x = 0; x < n; x++) 
	    	assignment[x]=xy[x];//forming answer there
	        //ret += cost[x][xy[x]];
	    return;
	}
	
	void update_labels()
	{
	    int x, y;
	    double delta = INF;             //init delta as infinity
	    for (y = 0; y < n; y++)            //calculate delta using slack
	        if (!T[y])
	            delta = Math.min(delta, slack[y]);
	    for (x = 0; x < n; x++)            //update X labels
	        if (S[x]) lx[x] -= delta;
	    for (y = 0; y < n; y++)            //update Y labels
	        if (T[y]) ly[y] += delta;
	    for (y = 0; y < n; y++)            //update slack array
	        if (!T[y])
	            slack[y] -= delta;
	}
	void add_to_tree(int x, int prevx, double[][] cost) 
	//x - current vertex,prevx - vertex from X before x in the alternating path,
	//so we add edges (prevx, xy[x]), (xy[x], x)
	{
	    S[x] = true;                    //add x to S
	    prev[x] = prevx;                //we need this when augmenting
	    for (int y = 0; y < n; y++)    //update slacks, because we add new vertex to S
	        if (lx[x] + ly[y] - (-cost[x][y]) < slack[y])
	        {
	            slack[y] = lx[x] + ly[y] - (-cost[x][y]);
	            slackx[y] = x;
	        }
	}
	
	void augment (double[][] cost)                         //main function of the algorithm
	{
		try{
	    if (max_match == n) return;        //check wether matching is already perfect
	    int x, y, root=-1;  
	    int[] q = new int[n];//just counters and root vertex
	    int wr = 0, rd = 0;          //q - queue for bfs, wr,rd - write and read
	                                       //pos in queue
	    Arrays.fill(S, false);       //init set S
	    Arrays.fill(T, false);          //init set T
	    Arrays.fill(prev, -1);
	    
	    for (x = 0; x < n; x++)            //finding root of the tree
	        if (xy[x] == -1)
	        {
	            q[wr++] = x;
	            root = x;
	            prev[x] = -2;
	            S[x] = true;
	            break;
	        }

	    for (y = 0; y < n; y++)            //initializing slack array
	    {
	        slack[y] = lx[root] + ly[y] - (-cost[root][y]);
	        slackx[y] = root;
	    }
	    while (true)                                                        //main cycle
	    {
	        while (rd < wr)                                                 //building tree with bfs cycle
	        {
	            x = q[rd++];                                                //current vertex from X part
	            for (y = 0; y < n; y++)                                     //iterate through all edges in equality graph
	                if ((-cost[x][y]) == lx[x] + ly[y] &&  !T[y])
	                {
	                    if (yx[y] == -1) break;                             //an exposed vertex in Y found, so
	                                                                        //augmenting path exists!
	                    T[y] = true;                                        //else just add y to T,
	                    q[wr++] = yx[y];                                    //add vertex yx[y], which is matched
	                                                                        //with y, to the queue
	                    add_to_tree(yx[y], x, cost);                              //add edges (x,y) and (y,yx[y]) to the tree
	                }
	            if (y < n) break;                                           //augmenting path found!
	        }
	        if (y < n) break;                                               //augmenting path found!

	        update_labels();                                                //augmenting path not found, so improve labeling
	        wr = rd = 0;                
	        for (y = 0; y < n; y++)        
	        //in this cycle we add edges that were added to the equality graph as a
	        //result of improving the labeling, we add edge (slackx[y], y) to the tree if
	        //and only if !T[y] &&  slack[y] == 0, also with this edge we add another one
	        //(y, yx[y]) or augment the matching, if y was exposed
	            if (!T[y] &&  slack[y] == 0)
	            {
	                if (yx[y] == -1)                                        //exposed vertex in Y found - augmenting path exists!
	                {
	                    x = slackx[y];
	                    break;
	                }
	                else
	                {
	                    T[y] = true;                                        //else just add y to T,
	                    if (!S[yx[y]])    
	                    {
	                        q[wr++] = yx[y];                                //add vertex yx[y], which is matched with
	                                                                        //y, to the queue
	                        add_to_tree(yx[y], slackx[y], cost);                  //and add edges (x,y) and (y,
	                                                                        //yx[y]) to the tree
	                    }
	                }
	            }
	        if (y < n) break;                                               //augmenting path found!
	    }

	    if (y < n)                                                          //we found augmenting path!
	    {
	        max_match++;                                                    //increment matching
	        //in this cycle we inverse edges along augmenting path
	        for (int cx = x, cy = y, ty; cx != -2; cx = prev[cx], cy = ty)
	        {
	            ty = xy[cx];
	            yx[cy] = cx;
	            xy[cx] = cy;
	        }
	        augment(cost);                                                      //recall function, go to step 1 of the algorithm
	    }
		}catch (Exception e){
			System.err.println("Exception: " + e.getMessage());
			e.printStackTrace();
		}
	}//end of augment() function 
	

}
