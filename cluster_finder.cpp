/* NOTA: Nei commenti a questo codice scriverò a fianco di cicli e invocazioni ricorsive non banali la complessità computazionale dell'esecuzione di tali cicli,
intendendo con n=#V ed m=#E il numero di nodi ed archi del grafo (ricordiamo che comunque per grafi reali #V~#E) e con #friends il numero di "amici" di un nodo,
esso è limitato e può essere ignorato ai fini di un'analisi asintotica, ai fini di un'analisi reale si può invece notare, dividendo m per n, che la sua media è
bassa (e non sarà molto diverso per gli altri momenti: E[#friends^2], E[#friends^3], etc...) dunque possiamo ignorare questi fattori nella nostra analisi della
complessità dell'algoritmo.

Infine otteniamo la complessità di m*log(m)*long(n) */

#include<iostream>
#include<fstream>
#include<vector>
#include<iomanip>
#include<algorithm>

//PARAMETRI DA CUI DIPENDE L'ALGORITMO

#define REDUCTION_FACTOR 0.1		// frazione di cluster accorpati ad ogni iterazione di Kruskal troncato.
#define THRESHOLD_FACTOR 0.9			// minima frazione dei dei nodi nei cluster da stampare.
#define MIN_PRINT_TO_ITERATE 5		// minimo numero di cluster stampabili affinchè si continui ad iterare Kruskal troncato.
#define PRINT_IF_LESS_THAN 20		// massimo numero di cluster adatti alla stampa (ovvero t.c. size>=THRESHOLD_FACTOR) affinchè vengano stampati su "modularity".
using namespace std;

typedef pair<int, int> iPair;

struct Node{						// struttura di cui utilizzeremo un vector per memorizzare il grafo, il flag ci dice se il nodo è ancora 
	int label;						// presente nel grafo durante la fase di preprocessing in cui suddividiamo in cliques il grafo.
	vector<int> friends;
	bool flag;
	Node(int l){label=l;}
	Node(){}
};
struct KruskalGraph{						// struttura che utilizzeremo per memorizzare il grafo quando vorremo applicarci Kruskal in quanto
	int V, E;								// ci è utile raccorgliere gli edges in un solo array per ordinarli con la funzione
	vector<pair<double, iPair>> edges;		// sort della STL.
	KruskalGraph(int x, int y){V=x; E=y;}
	KruskalGraph(){}
	void addE(long double d, int i, int j){
		edges.push_back(make_pair(d, make_pair(i,j)));
	}
};

KruskalGraph setForKruskal(vector<vector<int>> clusters, vector<Node> refer, vector<int> cliBelong);
bool heavier(pair<double, iPair> x, pair<double, iPair> y);
int find(int x);
void unite(int x, int y);
vector<vector<int>> clusterize(vector<vector<int>> cliques, vector<int>& cliBelong);
vector<vector<int>> filter(vector<vector<int>> clusters, int nodeN);
bool biggerFirst(vector<int> x, vector<int> y);
void printResult(vector<vector<int>> clusters, vector<Node> refer, int m);

vector<int> father;							// variabile globale utilizzata per denotare le componenti connesse in MST
ofstream file_statistics, file_belongings;	

int main(){	
	int i, j, nodeN=0, edgeN=0;;									//leggo da file e salvo in vector<Node> stone
		
	ifstream myFile;
	string fileName;
	cout<<"Inserire il nome del file"<<endl;
	cin>>fileName;
	myFile.open(fileName, ios::in);
	file_statistics.open("modularity", ios::out);
	file_belongings.open("clusters_belongings", ios::out);
	while(myFile>>i) if(i>nodeN) nodeN=i;
	nodeN++;
	myFile.clear();
	myFile.seekg(0, ios::beg);
	
	vector<int> dinamicTab(nodeN, -1), staticTab(nodeN, -1);
	vector<bool> stoneTab(nodeN, false);
	vector<Node> stone(nodeN);

	for(i=0; i<nodeN; i++){
		stone[i].label=i;
		stone[i].flag=true;
	}
	while(myFile>>i>>j){
		if(i!=j){								// NOTA: i!=j, ovvero non ammettere self-loops è un preprocessing necessario per non 
			stone[i].friends.push_back(j);		// mandare in loop l'algoritmo.
			stone[j].friends.push_back(i);
			edgeN++;			
		}
	}
									// prima parte: suddivido in cliques il grafo di partenza
	int survived=nodeN, pivot=0;	// survived conta i nodi residui nel grafo per il ciclo maggiore, pivot il più piccolo
	vector<vector<int>> cliques;	// dei nodi superstiti, si può eventualmente trovare un'euristica per la scelta di tale pivot,
	vector<Node> refer=stone;		// ho provato a scegliere il nodo con grado minore nel grafo residuo (idea presa dalla letteratura
	vector<int> cliBelong(nodeN);	// per ottenere una pseudo partizione minimale in cliques dei vertici) ma questo avrebbe aumentato ad
									// O(n^2) la complessità di questa fase senza produrre netti miglioramenti sulla qualità dei clusters finali.
									
	while(survived>0){ // O(n*#friends^3)=O(n)
		vector<int> cli;
		while(!stone[pivot].flag) pivot++;  // scelgo il minimo nodo ancora vivo
		cli.push_back(pivot);				// lo aggiungo alla mia clique
		
		vector<Node> alive;
		vector<int>::iterator ifr;
																	
		for(ifr=stone[pivot].friends.begin(); ifr!=stone[pivot].friends.end(); ifr++){ 		// O(#friends)
			if(stone[*ifr].flag){				// copio in alive gli amici del mio pivot
				alive.push_back(Node(*ifr));
				stoneTab[*ifr]=true;			// flaggo in stoneTab i nodi attivi nel mio grafo
			}
		}
		for(int i=0; i<alive.size(); i++){		//O(#friends^2)
			vector<int> tmp=stone[alive[i].label].friends; 
			for(int j=0; j<tmp.size(); j++){				//O(#friends)	
				if(stoneTab[tmp[j]]) alive[i].friends.push_back(tmp[j]); // inizializzo gli amici in alive
			}
		}
		for(vector<Node>::iterator it=alive.begin(); it!=alive.end(); it++){ // O(#friends)
			stoneTab[it->label]=false; // pulisco la mia struttura per segnare i nodi in alive in modo efficiente
		}
		// Il seguente loop itera ciò che ho appena fatto tra stone ed alive senza però utilizzare la struttura rigida di stone
		// per cui vale per ogni i stone[i].label==i, ma utilizzando un coppia di tabelle contenenti gli indirizzi (sotto forma
		// di indici dei vari elementi di alive.
		while(!alive.empty()){ // O(#friends^3) poichè la dimensione di alive è limitata da #friends
			cli.push_back(alive[0].label); 
			
			vector<Node> newAlive;
			
			for(int i=0; i<alive.size(); i++) staticTab[alive[i].label]=i; // O(#friends), scrivo gli indirizzi dei nodi in alive
			for(ifr=alive[0].friends.begin(); ifr!=alive[0].friends.end(); ifr++){ // O(#firends)
				dinamicTab[*ifr]=staticTab[*ifr];	// scrivo gli indirri in alive dei nodi in newAlive
				newAlive.push_back(Node(*ifr));
			}			
			for(int i=0; i<newAlive.size(); i++){ //O(#friends^2)
				vector<int> tmp=alive[dinamicTab[newAlive[i].label]].friends;
				for(int j=0; j<tmp.size(); j++){ // O(#friends)
					if(dinamicTab[tmp[j]]>=0) newAlive[i].friends.push_back(tmp[j]); // inserisco in newAlive gli amici
				}
			}			
			for(int i=0; i<alive.size(); i++) staticTab[alive[i].label]=-1; //, O(#friends), pulisco le mie tabelle in modo efficiente
			for(int i=0; i<newAlive.size(); i++) dinamicTab[newAlive[i].label]=-1;
			alive=newAlive;
		}
		for(ifr=cli.begin(); ifr!=cli.end(); ifr++){ // O(1), poichè la dimensione delle cliques è limitata (alternativamente si può
													// tra tutte le iterazioni si entra al più n volte in questo ciclo
			cliBelong[*ifr]=cliques.size(); // inizializzo cliBelong, vettore che mi permette di accedere, dato il label di un
			stone[*ifr].flag=false;			// nodo alla posizione della sua clique di appartenenza in tempo costante.
			survived--;
		}
		cliques.push_back(cli);		
	}
	// A questo punto il mio vettore cliques contiene dei vettori di labels rappresentanti la partizione in cliques del nostro grafo
	//Seconda parte: algoritmo di Kruskal troncato ed iterato
	
	vector<vector<int>> clusters=cliques;					
	KruskalGraph cnet;
	int iterNumber=0;	// conta le iterazioni eseguite
	int printable=MIN_PRINT_TO_ITERATE+1; // inizializzazione utile a verificare la prima guardia
	
	while(printable>MIN_PRINT_TO_ITERATE){ // O(m*log(m)*log(n) in quanto il ciclo viene eseguito log(n) volte: ad ogni passo riduco 
											// di una frazione costante il numero di nodi, itero finchè ho abbastanza nodi da stampare

		cnet=setForKruskal(clusters, refer, cliBelong); // data la suddivisione in clusters costruisce il grafo avente i clusters come nodi
		father.resize(cnet.V);							// e la frequenza dei collegamenti tra loro come archi pesati
		for(int i=0; i<cnet.V; i++) father[i]=i;		
		sort(cnet.edges.begin(), cnet.edges.end(), heavier); //O(m*log(m)), ordino gli archi dal più pesante al più leggero
		
		int mstEdges=0, mstCurr=0;
		
		while(mstEdges<cnet.V*REDUCTION_FACTOR && mstCurr<cnet.E){//O(n), itero fino ad accorpare la frazione di cluster REDUCTION_FACTOR
			int a=cnet.edges[mstCurr].second.first;	
			int b=cnet.edges[mstCurr].second.second;
			if(find(a)!=find(b)){ // se non appartengono alla stessa componente connessa (espressa da father),
				unite(a, b);	// aggiungo il nodo che li congiunge
				mstEdges++;
			}
			mstCurr++;
		}
		int oldSize=clusters.size();
		clusters= clusterize(clusters, cliBelong);	// accorpa i cluster appartenenti alle stesse componenti connesse di MST
		if(oldSize==clusters.size())  break; 		//interrompe il ciclo se il numero di cluster si stabilizza
		iterNumber++;
		cout<<"Iterazione "<<iterNumber<<endl;		// stampa su std output il numero dell'iterazione correntemente eseguita
		vector<vector<int>> tmp=clusters;
		tmp=filter(clusters, nodeN); // O(n*log(n)), inserisce solo i primi cluster,che insieme constano del THRESHOLD_FACTOR dei nodi
		if(tmp.size()<=PRINT_IF_LESS_THAN){ // se sono pochi abbastanza li stampo sui file di output
			// stampo la tabella sul file "modularity"
			file_statistics<<"Iterazione "<<iterNumber<<endl;
			// stampo sul file clusters_belongings per ogni nodo la coppia label-indice del cluster (questi sono ordinati per size decrescente)
			file_belongings<<"Iterazione "<<iterNumber<<endl;
			printResult(tmp, refer, edgeN);
		}
		printable=tmp.size();
		cout<<100*THRESHOLD_FACTOR<<"%"<<" dei nodi è contenuto nei primi "<<printable<<" clusters"<<endl<<endl;
	}	
	
	myFile.close();
	file_statistics.close();
	file_belongings.close();
	return 0;
}

bool biggerFirst(vector<int> x, vector<int> y){return x.size()>y.size();}

KruskalGraph setForKruskal(vector<vector<int>> clusters, vector<Node> refer, vector<int> cliBelong){ //O(n*#friends)
	KruskalGraph cnet(clusters.size(), 0);
	vector<int>::iterator in, it;
	
	for(int i=0; i<clusters.size(); i++){ // O(n*#friends)
		vector<int> weights(clusters.size(), 0);
		for(in=clusters[i].begin(); in!=clusters[i].end(); in++){
			for(it=refer[*in].friends.begin(); it!=refer[*in].friends.end(); it++){
				weights[cliBelong[*it]]++;
			}
		}
		for(int j=i+1; j<clusters.size(); j++){ //O(n)
			if(weights[j]!=0){
				cnet.addE((double) weights[j]/(2*clusters[i].size()*clusters[j].size()), i, j);
				cnet.E++;
			}
		}
	}	
	return cnet;
}

bool heavier(pair<double, iPair> x, pair<double, iPair> y){
	return x.first>y.first;
}

void unite(int x, int y){ // O(n), unisce le componenti connesse di x ed y
	int fx=find(x);
	int fy=find(y);
	father[fx]=fy;
}

int find(int x){ // O(n), trova il padre di x, il quale identifica la componente connessa di ogni elemento
	if(father[x]==x) return x;
	return find(father[x]);
}

vector<vector<int>> clusterize(vector<vector<int>> cliques, vector<int>& cliBelong){ // O(n)
	vector<vector<int>> clusters;
	vector<int> map(cliques.size(), -1);
	vector<int>::iterator in;
	
	for(int j=0; j<cliques.size(); j++){ //O(n)
		if(map[find(j)]<0){
			map[find(j)]=clusters.size();
			clusters.push_back(vector<int>());
		}
		for(in=cliques[j].begin(); in!=cliques[j].end(); in++){ //O(n)
			clusters[map[find(j)]].push_back(*in);
			cliBelong[*in]=map[find(j)];
		}
	}
	return clusters;
}

vector<vector<int>> filter(vector<vector<int>> clusters, int nodeN){ // O(n*log(n))
	vector<vector<int>>::iterator it;	// seleziona, una volta ordinati, i primi cluster per grandezza finchè
	vector<vector<int>> cl2;			// la somma delle loro carinalità non eccede la frazione soglia dei nodi totali
	int acc=0;
	sort(clusters.begin(), clusters.end(), biggerFirst); // O(n*log(n))
	for(it=clusters.begin(); it!=clusters.end() && (acc+=it->size())<nodeN*THRESHOLD_FACTOR; it++){ // O(n)
		cl2.push_back(*it);
	}
	return cl2;
}

void printResult(vector<vector<int>> clusters, vector<Node> refer, int m){ //ATTENZIONE! Il file sarà leggibile impostando 1 tab = 8 spazi.
	double weights[clusters.size()][clusters.size()];					// Tale impostazione è stata scelta in quanto in fase di test era quella
	vector<int> cliBelong(refer.size(), -1);							// visualizzata correttamente da terminale, inoltre è l'impostazione di 
	vector<Node>::iterator in;											// default del mio lettore, sperabilmente anche di altri. 
	vector<int>::iterator ifr;
	int degSum[clusters.size()];
	
	//NOTA: questa funzione viene invocata solo quando clusters.size()<PRINT_IF_LESS_THAN, dunque tutti i cicli su tali indici sono trascurabili
	// In totale questa funzione ha complessità O(n*#friends)
	
	for(int i=0; i<clusters.size(); i++){
		degSum[i]=0;
		for(int j=0; j<clusters.size(); j++){
			weights[i][j]=0;
		}
	}
	for(int i=0; i<clusters.size(); i++){ // O(n)
		for(ifr=clusters[i].begin(); ifr!=clusters[i].end(); ifr++){ // O(n)
			cliBelong[*ifr]=i;
		}
	}	
	for(in=refer.begin(); in!=refer.end(); in++){ // O(n*#friends)
		for(ifr=in->friends.begin(); ifr!=in->friends.end(); ifr++){ //O(#friends)
			if(cliBelong[in->label]>=0 && cliBelong[*ifr]>=0) weights[cliBelong[in->label]][cliBelong[*ifr]]++;
			if(cliBelong[in->label]>=0) degSum[cliBelong[in->label]]++;
			if(cliBelong[*ifr]>=0) degSum[cliBelong[*ifr]]++;
		}
	}
	//divido per due gli elementi di degSum e della diagonale di weights poichè contati 2 volte
	for(int i=0; i<clusters.size(); i++){ 
		weights[i][i]=weights[i][i]/2;
		degSum[i]=degSum[i]/2;
	}
	// stampo la tabella sul file "modularity"
	file_statistics<<endl;
	for(int i=0; i<clusters.size(); i++) file_statistics<<"\t\t"<<clusters[i].size(); 
	file_statistics<<endl<<endl;
	for(int i=0; i<clusters.size(); i++){
		file_statistics<<clusters[i].size();
		for(int j=0; j<clusters.size(); j++){
			if(weights[i][j]!=0){
				double tmp=(double) weights[i][j]/m;
				double tmp2=(double)degSum[i]/(2*m)*(double)degSum[j]/(2*m);
				tmp-=tmp2;
				file_statistics<<"\t"<<setprecision(8)<<tmp;
			}
			else file_statistics<<"\t"<<0<<"\t";
		}
		file_statistics<<endl;
	}
	file_statistics<<endl<<endl;
	// stampo sul file clusters_belongings per ogni nodo la coppia label-indice del cluster (questi sono ordinati per size decrescente)
	// NOTA: sono stampati solo i nodi appartenenti ai clusters visualizzati in "modularity"
	for(vector<Node>::iterator it=refer.begin(); it!=refer.end(); it++){ // O(n);
					if(cliBelong[it->label]>=0) file_belongings<<it->label<<" "<<cliBelong[it->label]<<endl;
				}
				file_belongings<<endl;
}
