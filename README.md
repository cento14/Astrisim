# Astrisim

Simulatore scientifico del MiniArray di ASTRI
Installazione :

Scompattando il file .tar viene creata una directory astrisim/ con il codice e 3 sottodirectory :

responses/, che contiene i file che descrivono le risposte dello strumento, 
sources/, che contiene i file che descrivono le sorgenti e 
simulations/, che contiene l'output del simulatore

AstriSim richiede che siano installate le astrolib di IDL.

Istruzioni per l'utilizzo :

Far partire IDL nella directory astrisim/ e dare il comando: 

IDL> astrisim 

(nota: se si vuole usare una directory diversa occorre aggiungere a !path la directory scelta) 

Gli input, output e parametri vari possono essere modificati attraverso il parameter file : astrisim.par. In particolare :

ins_path e ins_file : descrivono la directory e il nome del file che descrive la risposta dello strumento. Il formato di questi file e' uguale a quello di ctobssim , il simulatore di CTA. (Il file usato di default e' array7w3.dat, con i valori prodotti da Stamerra et al.)
src_path e src_list : sono la directory e la lista dei file che descrivono le sorgenti che si vogliono includere nella simulazione (a titolo di esempio sono gia' presenti i file relativi a Crab e SNR IC443)
t_hours : durata dell'osservazione in ore
l_asse e b_asse : direzione del centro del campo di vista in coordinate galattiche
out_path e out_file : directory e nome del file di output
Output :

L'output della simulazione e' un file fits che contiente la events-list (lista di direzione ricostruita, energia ricostruita e tempo) degli eventi dati dalle sorgenti e background. Il formato del file e' compatibile con i CTOOLS. Le osservazioni simulate con AstriSim possono quindi essere analizzate/visualizzate/manipolate con i CTOOLS o altri software comunemente usati (ds9, fv, ...). Di particolare utilita' sono ctbin, che permette di creare mappe di cielo in diversi range di energia e ctlike che fa un'analisi likelihood sulle sorgenti (trova flussi, spettri, posizioni...).

Sorgenti :

Le sorgenti sono descritte dai file nella directory sources, ogni file descrive una sorgente. Nel file sono descritte le coordinate celesti della sorgente e il flusso nelle varie bande di energia. In futuro sara' possibile, nello stesso file includere anche una curva di luce. I file sono in formato fits, editabili con fv o con sw simili. 
Questi file possono essere crati anche con la routine mksrc inclusa in AstriSim. Per esempio con il comando : 

IDL> mksrc, 'pippo.fits', 100, -3, index=-2, k=1e-11 

viene prodotto il file pippo.fits che descrive una sorgente posta in l,b = 100., 3 e con spettro F = k * (E / 1 TeV ) ^(index) espresso in : fotoni / cm^2 sec TeV.

Un esempio :

Lanciando il simulatore senza modificare il file astrisim.par gia' presente in astrisim.tar.gz, viene simulato un'osservazione dell'anticentro Galattico della durata di 30 ore, con le sorgenti Crab e IC443 nelle rispettive posizioni. Si ottiene: 

______________ ASTRI SCIENTIFIC SIMULATOR (v2.0) _____________ 

Instrument Response File : responses/array7w3.dat 
Simulated Sources : sources/crab.fits sources/ic443.fits 

Livetime (hours) : 30 
Pointing Direction (l,b) : 187.000 -2.00000 

Sources photons (plotted in red) : 55311 
Total Background events : 25462 

Saving file... (simulations/anticentro.evt) 
Done! 

e un plot che da la distribuzione in cielo degli eventi prodotti (quelli delle sorgenti sono in rosso). Nella directory simulations/ viene salvato il file anticentro.evt che contiene ~80000 eventi. 
Usando ctbin si possono quindi produrre le mappe di cielo mostrate sotto, e usando ctlike analizzare le due sorgenti. (NB. le sorgenti sono, in questo esempio, modellizzate da delle semplici power-law. Crab ha in realta' un cutoff intorno a 15 TeV. Il flusso sopra 15 TeV e' quindi sovrastimato. ) 

