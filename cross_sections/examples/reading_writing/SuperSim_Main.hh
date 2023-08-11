#ifndef SuperSim_Main_hh
#define SuperSim_Main_hh 1
// $Id: SuperSim_Main.hh,v 1.10 2015/07/07 20:21:18 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//  File:        SuperSim_Main.hh                                     //
//  Description: Driver for SuperCDMS GEANT4 simulation.  If physics  //
//		 other than the default is needed, must use default   //
//		 constructor and call Configure(), Run() explicitly.  //
//                                                                    //
//  Author:      Michael Kelsey (SLAC)                                //
//  Date:        09 May 2012                                          //
//                                                                    //
//  20120712  Add public SetVerboseLevel() to pass through to user    //
//	      actions, which will be accessed via RunManager.         //
//  20120718  Add SuperSim_Messenger to get verbose pre-defined.      //
//  20130523  Remove SetRunNumber function (users must migrate to new //
//            macro command); add function for command-line aliases.  //
//  20140613  Use custom RunManager to defer execution of RM->Init()  //
//  20150706  Replace direct registration of physics list with use of //
//		CDMSPhysListFactory.                                  //
//  20151221  For G4 10.x, move action registration to new class      //
//  20200310  Add static function to print out version on request.    //
//  20200511  Follow move of CDMSRunManager to CDMSactions.           //
//  20230810  A. Biffl -- this version of the header file is being    //
// 		used in this project to allow access to previously    //
// 		private member variables                              //
////////////////////////////////////////////////////////////////////////

#include "globals.hh"
#include "CDMSactions/CDMSRunManager.hh"

class CDMSActionInitialization;
class CDMSGeomConstructor;
class G4UImanager;
class G4VUserPhysicsList;
class G4VisManager;
class SuperSim_Messenger;


class SuperSim_Main {
public:
  SuperSim_Main();		// No actions, user must Configure(), Run()
  SuperSim_Main(int argc, char* argv[]);	// Use default physics
  ~SuperSim_Main();

  void Configure(const G4String& plName="");	// Specify custom physics list
  void Run(int argc, char* argv[]);		// Pass command-line arguments
  void Run(const char* macro=0);		// Specify macro file directly

  void SetVerboseLevel(G4int verbose);		// Pass verbosity to actions

  // Special function to print version info and exit, if requested
  static void VersionCheck(int argc, char* argv[]);

  void SetDirectory();			// Set macro alias pointing to base
  void SetArguments(int argc, char* argv[]);	// Copy arguments to aliases

  CDMSRunManager runManager;		// This must be instantiated early!
  G4VisManager* visManager;		// Will be instantiated in Configure
  G4UImanager* UImanager;		// Singleton, not owned by SuperSim

  // Cache user actions for verbosity management
  // NOTE:  These objects are NOT owned, RunManager will delete them
  CDMSActionInitialization* act;
  CDMSGeomConstructor* geo;

  SuperSim_Messenger* messenger;
};

#endif	/* SuperSim_Main_hh */
