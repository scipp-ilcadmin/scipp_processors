#undef _GLIBCXX_USE_CXX11_ABI
#define _GLIBCXX_USE_CXX11_ABI 0
#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <random>

#include "lcio.h"
#include "IMPL/LCEventImpl.h"
#include "EVENT/LCEvent.h"
#include "EVENT/LCCollection.h"
#include "EVENT/SimCalorimeterHit.h"
#include "UTIL/CellIDDecoder.h"
#include "IMPL/LCRunHeaderImpl.h"


#include "include/simple_list_geometry_edit.h"
#include "include/beamcal_scanner_edit.h"
#include "include/beamcal_reconstructor_edit.h"
#include "include/BeamCalRecon_edit.h"

#include "scipp_ilc_globals.h"


using namespace std;

namespace scipp_ilc {
    namespace beamcal_recon_C {


/*
 * In the lcdd file, the beamcal has a square pixel layout. But we want
 * a radial pixel layout. So we have the pixelate_beamcal function to 
 * change the pixelation. Unfortunately, all the hits on a pixel are 
 * reduced to a single point in the center of the pixel. So that 
 * function also acts to "distribute" the single hit evenly to several
 * hits across the square pixel area, essentially acting as an
 * anti-aliasing function. _spreadfactor denotes how heavily we are anti-
 * aliasing (how much we are dividing up a single hit). If _spreadfactor
 * = 1, then no spreading is performed. With 2, the hit is divided into
 * a 2x2 hits; with 3, you get 9 hits in a 3x3 pattern; etc.
 * _cellsize is the original square pixel size as defined in the lcdd and 
 * compact.xml.
 */

        static const float _cellsize = 1; //milimeter
        static const float _spreadfactor = 1; //1; we decided we don't need to spread a 1 mm pixel
        static const bool _remove_negative = true;

        //the fraction of background events that the program is
        //allowed to reject. This is used to calculate the
        //sigma cut
        static const float _rejection_limit = 0.1;

        static float _sigma_cut;

        static bool _adding_to_stats;
        static int _num_bgd_events;
        static unordered_map<int,double>* _energy_totals;
        static unordered_map<int,double>* _square_energy_totals;
        static unordered_map<int,double>* _energy_averages;
        static unordered_map<int,double>* _energy_std_devs;
        static unordered_map<int,int>* _times_hit;
        
        vector<pixel_map*>* _database;



        /*
         * Open an lcio event and take the rectilinear beamcal hits and apply them
         * to the radial tiling scheme we use. If the simulation-level pixel size
         * is too large, you also need to use the hit spreader.
         *
         * Additionally, if processing the background hits,
         * load up the statistics maps.
         *
         *
         * IMPORTANT: From this point forward, all "pixels" are NOT just one pixel.
         * They are "layer compressed" pixels. That is, all pixels of the same ID
         * between layers "layer_min" and "layer_max" are all compressed together
         * into a single pixel. This drastically reduces processing time and memory
         * consumption, and also aids identification of signal events, as signal
         * events penetrate more deeply into the beamcal than background events.
         * To understand why we choose to between layer 10 and layer 40 (here
         * layer 9 and 39, because indexing starts at layer 0), please refer to
         * Alex Bogert's Thesis.
         *
         * Developer Note: A possible improvement on the current algorithm may
         * be to weight hits from later layers. That is, a layer 9 hit simply
         * adds hit->getEnergy(), but a layer 28 hit will add:
         * layerWeight(28)*hit->getEnergy(), where layerWeight(int layer) is
         * a function you would need to define that fits the signal's energy
         * deposition per layer. The signal's energy deposition per layer is
         * higher than bgd events, so magnifying late-layer energy deposition
         * may increase signal recognition. Of course, this will make it more
         * of a pain to reconstruct the signal event's energy later on, once
         * you've identified it (you'll probably just have to rerun over the
         * beamcal hits and only accept the hits that land within the IDs
         * specified in the cluster IDlist returned by the scanner). 
         *
         */
        static void pixelate_beamcal(lcio::LCEvent* event, pixel_map* new_pixels) {
            double dim = _cellsize / ( _spreadfactor );
            double Ediv = (_spreadfactor * _spreadfactor);

            unsigned int layer_min = 9;
            unsigned int layer_max = 39;
	    

            lcio::LCCollection* col = event->getCollection("BeamCalHits") ;
            if( col != NULL ){
                lcio::CellIDDecoder<lcio::SimCalorimeterHit> decoder = lcio::CellIDDecoder<lcio::SimCalorimeterHit>(col);


                int nElements = col->getNumberOfElements()  ;
		float old_z;
		float old_y;
		float old_x;
		float radius;
		float old_energy;
		unsigned int layer;
		float spread_energy;
		float spread_x;
		float spread_y;
		int ID;
		float energy;

                for(int hitIndex = 0; hitIndex < nElements ; hitIndex++){
                    lcio::SimCalorimeterHit* hit = dynamic_cast<lcio::SimCalorimeterHit*>( col->getElementAt(hitIndex) );



                    const float* old_pos = hit->getPosition();

                    old_z = old_pos[2];
                    old_y = old_pos[1];
                    old_x = old_pos[0] - abs(old_z)*_transform;
                    radius = hypot(old_x,old_y);
                    old_energy = hit->getEnergy();
                    layer = decoder(hit)["layer"];

                    if ( _remove_negative && (old_z<0) ) continue;
                    if ( radius > _radius_cut ) continue;
                    if ( layer < layer_min or layer_max < layer ) continue;
                    
                    if (_spreadfactor > 1) { // _spreadfactor == 1 , therefore fails condition
                        spread_energy = old_energy / Ediv;
                        for (int i = 0; i < _spreadfactor; i++) {
                            spread_x = (i*dim) + old_x + (dim/2.0) - (_cellsize/2.0);
                            for (int j = 0; j < _spreadfactor; j++) {
                                spread_y = (j*dim) + old_y + (dim/2.0) - (_cellsize/2.0);
                                ID = getID(spread_x,spread_y);
                                (*new_pixels)[ID] += spread_energy;
                            }
                        }
                    } else {
                        ID = getID(old_x,old_y);
                        (*new_pixels)[ID] += old_energy;
			//			cout << " old x:"<< old_x << ", old y:"<< old_y << endl;
			if(_adding_to_stats==true){
			  //			  _hitmap_bgd->Fill(old_x,old_y);


			}
			else{
			}
                    }
                }

                if (_adding_to_stats) {
                    for( auto pixel : *new_pixels ) {
                        ID = pixel.first;
                        energy = pixel.second;

                        (*_energy_totals)[ID] += energy;
                        (*_square_energy_totals)[ID] += ( (double)energy ) * ( (double)energy );
                        (*_times_hit)[ID] += 1;
                    }
                }
            }
        }



        /*
         * Read in the background file list, iterate through it line
         * by line, read in each slcio file, read the slcio's event by
         * event, and load the beamcal hits from each event into a 
         * pixel_map for that specific event.
         */
        static void process_background_events(string bgd_list_file_name) {
            _database = new vector<pixel_map*>();

            int numEventsRead = 0;
            try { 
                //open filelist
                ifstream filelist (bgd_list_file_name, ifstream::in);
                lcio::LCReader* lcReader = lcio::LCFactory::getInstance()->createLCReader() ;
                string slcioFile;
                lcio::LCEvent* event = NULL;

                //for each slcio file in the file list
                while ( filelist >> slcioFile ) {
                    lcReader->open(slcioFile);
                    
                    //for each event in the slcio file
                    while( (event=lcReader->readNextEvent()) ) {
                        pixel_map* new_pixels = new pixel_map();
                        pixelate_beamcal(event,new_pixels);
                        _database->push_back(new_pixels);

                        numEventsRead++;
                        cout << "Database read number = " << numEventsRead << endl;
                        if ( numEventsRead >= _num_bgd_events ) {
                            delete event;
                            break;
                        }
                    }
                    lcReader->close();
                    if ( numEventsRead >= _num_bgd_events ) break;
                }
                filelist.close();
            } catch(lcio::IOException& e) {
                cout << " Unable to read and analyze the LCIO file - " << e.what() << endl ;
            }
        }



        /*
         * Read in all of the bgd events, store their beamcal hit
         * information in the _database vector, and get the averages
         * and standard deviations of all the pixels over all events.
         */
        static void generate_database(string bgd_list_file_name) {
            cout << "Generating Database...\n";

             _energy_totals = new unordered_map<int,double>();
             _square_energy_totals = new unordered_map<int,double>();
             _energy_averages = new unordered_map<int,double>();
             _energy_std_devs = new unordered_map<int,double>();
             _times_hit = new unordered_map<int,int>();

            //read in all of the background events in the given
            //bgd file list.
            //NOTE: Depending on the number of bgd events, this
            //one function will take longer than the entire rest
            //of the reconstruction.
            process_background_events(bgd_list_file_name);
	    int ID;
	    double energy_total;
	    int hitcount;
	    double squared_energy_total;
	    double energy_average;
	    double average_of_squares;
	    double square_of_averages;
	    double energy_std_dev;

            for ( auto pixel : *_energy_totals ) {
                ID = pixel.first;
                energy_total = pixel.second;

                hitcount = (*_times_hit)[ID];
                squared_energy_total = (*_square_energy_totals)[ID];

                energy_average = energy_total / _num_bgd_events;

                average_of_squares = squared_energy_total / _num_bgd_events;
                square_of_averages = energy_average * energy_average;
                energy_std_dev = sqrt(average_of_squares - square_of_averages);

                if (hitcount == 1) {
                    energy_std_dev = -1.0;
                }

                (*_energy_averages)[ID] = energy_average;
                (*_energy_std_devs)[ID] = energy_std_dev;
            }


            delete _energy_totals;
            delete _square_energy_totals;
            delete _times_hit;

            cout << "Database succesfully generated.\n";
        }


        
        /*
         * Sort from greatest to least significance.
         * Largest significance cluster is zeroth element in list
         */
        static bool compare_cluster( beamcal_cluster* first, beamcal_cluster* second ) {
            return ( first->significance > second->significance );
        }



        /*
         * Run the clustering/signal identification algorithm for every
         * bgd event stored in the _database. This gives us the highest
         * significance cluster for each bgd event. We take each of these
         * clusters, and order them by their significance. Finally, we use
         * this sorted list to determine the sigma cut which only a certain
         * fraction (given by _rejection_limit) of the clusters (and thus
         * of the bgd events themselves) will exceed.
         */
        static void calibrate_scanner() {
            cout << "Calibrating Scanner...\n";

            vector<beamcal_cluster*> cluster_list;
            int map_num = 0;
            for( pixel_map* map : *_database ) {
                beamcal_cluster* new_cluster;
                new_cluster = scan_beamcal(map,_energy_averages,_energy_std_devs);
                cluster_list.push_back(new_cluster);
                cout << "   Calibrating on background event " << map_num++ << endl;
            }

            sort(cluster_list.begin(), cluster_list.end(), compare_cluster);
            
            int cutoff_index = (int)( cluster_list.size()*_rejection_limit );
            _sigma_cut = cluster_list[cutoff_index]->significance;

            cout << "Scanner calibration complete.\n";
        }



        /*
         * This function does three things: 
         * > setup the geometry,
         * > read in all the background events and setup the base statistics,
         * > get the signal sigma cut (the significance a cluster must have
         *          in order to be called a signal event)
         */
        void initialize_beamcal_reconstructor(string geom_file_name, string bgd_list_file_name, int bgd_events_to_be_read) {
            _num_bgd_events = bgd_events_to_be_read;

            //The _adding_to_stats variable is turned on during this sequence only.
            //It is only used only in the pixelate_beamcal function. However, that
            //function is used here, and also when processing signal events, and
            //when processing signal events you don't want to screw with the statistics.
            _adding_to_stats = true;
            initialize_geometry(geom_file_name); //from simple_list_geometry.h
            generate_database(bgd_list_file_name);
            calibrate_scanner();
            _adding_to_stats = false;

        }



        /*
         * Attempt to identify the cluster which marks the location of the
         * signal in this event. This function takes in an event which contains
         * a signal, and then makes a copy of a background event from the 
         * _database. The signal event is overlayed on top of the bgd event,
         * and then the scanner is invoked.
         */
        beamcal_cluster* reconstruct_beamcal_event(lcio::LCEvent* signal_event) {
            //literally copy-pasted this RNG from stack exchange
            //no idea how it works, but it does the job
            std::random_device rd; // only used once to initialise (seed) engine
            std::mt19937 rng(rd()); // random-number engine used (Mersenne-Twister in this case)
            std::uniform_int_distribution<int> uni(0,_num_bgd_events-1); // guaranteed unbiased
            int bgd_index = uni(rng);

            pixel_map bgd_populated_beamcal = *( (*_database)[bgd_index] );
            pixelate_beamcal( signal_event, &bgd_populated_beamcal );
            beamcal_cluster* signal_cluster = scan_beamcal(&bgd_populated_beamcal,_energy_averages,_energy_std_devs);

            signal_cluster->exceeds_sigma_cut = signal_cluster->significance > _sigma_cut;
            return signal_cluster;
        }
    }
}
