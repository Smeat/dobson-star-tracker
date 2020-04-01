/*
 * Conversion.cpp
 *
 *  Created on: 03.08.2019
 *      Author: lukas
 */
#include <Arduino.h>
#include <AccelStepper.h>
#include <TimeLib.h>

#include "./config.h"
#include "./macros.h"
#include "./conversion.h"
#include "./Observer.h"
#include "./location.h"

#ifdef SERIAL_DISPLAY_ENABLED
	#include "./display_unit.h"
#endif

// These are the coordinates that the telescope initially thinks it's pointed at.
// TODO This needs to be replaced by. Maybe by a config variable or depending on the opmode?
const double ra_h = 0;// 16.0;
const double ra_m = 0;// 41.7;
double ra_deg = (ra_h + ra_m / 60.0) * 15.0;

const double dec_d = 0;// 36.0;
const double dec_m = 0;// 28.0;
double dec_deg = dec_d + dec_m / 60.0;

RaDecPosition futurePosition = { 0, 0 };

boolean isPositiveDeclination = false;

boolean isAligned = false;

char txAR[10]; // Gets reported to stellarium when it asks for right ascension. Example: "16:41:34#"
char txDEC[11]; // Same as above with declination. Example: "+36d%c28:%02d#"

const byte numChars = 32;
char receivedChars[numChars];

bool newData = false; // Gets set to true whenever a complete command is buffered
static boolean recvInProgress = false; // True while a command is being received
static byte ndx = 0; // Number of command character received
const char startMarker = ':'; // Commands begin with this character
const char endMarker = '#'; // Commands end with this character


// These are the positions that can be cycled through with the position switch
// TODO Depending on the selected MOUNT_TYPE these need to change
// TODO If a diplay unit is connected, maybe these could be transferred somehow? Need to save memory on the nano though
// TODO Load from SD card?
const double debugPositions[][2] = {
	{ 250.42, 36.46 },
	{ 240.42, 35.46 },
	{ 230.42, 34.46 },
	{ 220.42, 33.46 },
	{ 210.42, 32.46 },
    { 180.42, 25.46 },
    { 150.42, 20.46 },
    { 100.42, 15.46 },
	{ 19.7, 36.5 },
	{ 17.7, 38.5 },
	{ 17.7, 39.5 },
	{ 279.2354166, 38.78555556 }, // Vega
	{ 213.9083333 , 19.17038889 },    // Arktur
	{ 37.9624166 , 89.2642777 }    // Polaris
};

const char* positionNames[] = {
	"Start",
	"Polaris",
	"Vega",
	"Arktur",
	"Pos 4",
	"Pos 5",
	"Pos 6",
	"Pos 7",
	"Pos 8",
	"Pos 9"
};

const int maxDebugPos = sizeof(debugPositions) / (sizeof(debugPositions[0]));


// If there is a target select button we need some variables
#ifdef TARGET_SELECT_PIN
	// Whether the target select button is pressed (true) or released (false)
bool targetSelectButtonPressed = false;

// The currently selected target position from the debugPositions array
// If the value is -1 it either means that no target is selected
// or that the target was selected differently (e.g. via :Sr...# and :Sd...# commands)
int selectedDebugTargetIndex = -1;

// This value tracks the last time the the target was changed.
// It is needed because we want to ignore repeated button presses
unsigned long timeOfLastTargetChange = 0;
#endif

void getRightAscension(Mount& scope);
void getDeclination(Mount& scope);

// This gets called by the Arduino setup() function and sends the initial position to Stellarium
void initCommunication(Mount& telescope) {
	getRightAscension(telescope);
	getDeclination(telescope);
}

/*
 * Every loop iteration this function checks if a new character is available from the Serial interface. It only reads the
 * character if the variable newData is set to false. It then tries to read a complete command (1 char per loop iteration).
 * When a command is complete, this function sets newData to true. Later on in the loop iteration the parseCommands function,
 * which tries to parse and execute the command, gets called. After command execution the newData variable is reset to false
 */
void receiveCommandChar() {
	if (Serial.available() > 0 && newData == false) {
		// Read the next character
		char rc = Serial.read();

		// Here we check if the command start marker has already been read in a previous loop iteration.
		if (recvInProgress == true) {
			// If the character is NOT the end marker, we treat it as part of the command, regardless of what it is
			if (rc != endMarker) {
				// Here we store the next character in the receivedChars array.
				receivedChars[ndx] = rc;
				ndx++;
				// This prevents the receivedChars array from overflowing
				if (ndx >= numChars) {
					ndx = numChars - 1;
				}
			} else {
				// End marker received, so the command is completely stored in receivedChars
				// Set the last value to \0 to terminate the string
				receivedChars[ndx] = '\0';
				// Set recvInProgress to false, so that next time a character is received we check for the start marker again
				recvInProgress = false;
				// Reset ndx to 0 so that next time we read a command it gets stored at the start of the receivedChars array again
				ndx = 0;
				// This signals to the parseCommands function that a complete command is now stored in receivedChars.
				newData = true;
			}
		} else if (rc == startMarker) {
			// Start marker received. Set recvInProgress to false, so that next time we check for the end marker or a command character
			recvInProgress = true;
		}
	}
}



// These will eventually be placeholders for a better looking command parser. Currently these do nothing
#define SERIAL_CMD_INVALID_COMMAND     -1
#define SERIAL_CMD_GET_RIGHT_ASCENSION 0
#define SERIAL_CMD_GET_DECLINATION     1
#define SERIAL_CMD_MOVE_STOP           2
#define SERIAL_CMD_MOVE_START          3
#define SERIAL_CMD_SET_RIGHT_ASCENSION 4
#define SERIAL_CMD_SET_DECLINATION     5




/**
 * Serial commands start
 * This section contains a function for each serial command that can get handled
 * TODO Refactor this into a class so that we can later support connections to different programs / tools
 */
// Print the possible commands
void printHelp() {
	Serial.println(":HLP# Print available Commands");
	Serial.println(":GR# Get Right Ascension");
	Serial.println(":GD# Get Declination");
	Serial.println(":Sr,HH:MM:SS# Set Right Ascension; Example: :Sr,12:34:56#");
	Serial.println(":Sd,[+/-]DD:MM:SS# Set Declination (DD is degrees) Example: :Sd,+12:34:56#");
	Serial.println(":MS# Start Move; Starts tracking mode if not enabled");
	Serial.println(":TRK0# Disable tracking");
	Serial.println(":TRK1# Enable tracking");
	Serial.print(":DBGM[0-" + String(maxDebugPos - 1) + "]# Move to debug position X");
	Serial.println(":DBGMIA# Increase Ascension by 1 degree");
	Serial.println(":DBGMDA# Decrease Ascension by 1 degree");
	Serial.println(":DBGMID# Increase Declination by 1 degree");
	Serial.println(":DBGMDD# Decrease Declination by 1 degree");
	Serial.println(":DBGDM[00-99]# Disable Motors for XX seconds");
	Serial.println(":DBGDSP# Send status update to display / serial console");
}



// Reports the current right ascension
void getRightAscension(Mount& scope) {
	const double pos = scope.getCurrentPosition().rightAscension / 15;

	const int hrs = static_cast<int>(pos);
	const int mins = static_cast<int>((pos - hrs) * 60);
	const int secs = static_cast<int>((pos - (hrs + mins / 60.)) * 3600);

	sprintf(txAR, "%02d:%02d:%02d#", hrs, mins, secs);

	Serial.print(txAR);
}

// Reports the current declination
void getDeclination(Mount &scope) {
	const double pos = scope.getCurrentPosition().declination;

	const int deg = static_cast<int>(pos);
	const int mins = static_cast<int>((pos - deg) * 60);
	const int secs = static_cast<int>((pos - (deg + mins / 60.)) * 3600);

	sprintf(txDEC, "%c%02d%c%02d:%02d#", deg > 0 ? '+' : '-', deg, 223, mins, secs);

	Serial.print(txDEC);
}

// Quit the current move by setting the target to the current position.
// This does not enable/disable tracking (if supported)
void moveQuit(Mount& scope) {
	// This sets the target to the position the telescope is actually pointing at, thus stopping any currently running moves
	scope.setTarget(scope.getCurrentPosition());
}

// Start the requested move
bool moveStart(Mount& scope) {
	// Immediately confirm to Stellarium
	Serial.print("0");

	// TODO Homing code needs to be better. It has to disable the steppers and there must be some way to enable/disable it
	// If homing mode is true we set isAligned to true
	// and return true to indicate to the loop() function that homing is complete.
	if (scope.getMode() == Mode::ALIGNING) {

		DEBUG_PRINTLN("Setting alignment");
		scope.setAlignment(futurePosition);
		isAligned = true;
		return true;
	}
	else {
		DEBUG_PRINTLN("Starting Move");
		scope.setTarget(futurePosition);
	}

	return false;
}


// Set Right Ascension (in hours, minutes and seconds)
// This doesn't yet set it on the telescope (happens in moveStart())
void setRightAscension(Mount &telescope) {
	// Immediately confirm to Stellarium
	Serial.print("1");

	// Parse the coordinates part of the command to integers
	const int hrs = multi_char_to_int(receivedChars[3], receivedChars[4]);
	const int mins = multi_char_to_int(receivedChars[6], receivedChars[7]);
	const int secs = multi_char_to_int(receivedChars[9], receivedChars[10]);

	DEBUG_PRINTLN();
	DEBUG_PRINTLN("Changing RA");
	DEBUG_PRINTLN(ra_deg);
	ra_deg = (hrs + mins / 60. + secs / 3600.) * 15;
	DEBUG_PRINTLN(ra_deg);

	// Store the new target right ascension
	futurePosition.rightAscension = ra_deg;
	
	// If there is a target select pin we need to reset the selected position to "none"
	#ifdef TARGET_SELECT_PIN
		selectedDebugTargetIndex = -1;
	#endif
}

// Set target Declination (in +/- degrees, minutes and seconds)
// This doesn't yet set it on the telescope (happens in moveStart())
void setDeclination(Mount &telescope) {
	// Immediately confirm to Stellarium
	Serial.print("1");

	// Whether the coordinates are positive (1) or negative (-1)
	const int multi = (receivedChars[3] == '+') ? 1 : -1;

	// Parse the coordinates part of the command to integers
	const int deg = multi_char_to_int(receivedChars[4], receivedChars[5]);
	const int mins = multi_char_to_int(receivedChars[7], receivedChars[8]);
	const int secs = multi_char_to_int(receivedChars[10], receivedChars[11]);

	DEBUG_PRINTLN();
	DEBUG_PRINTLN("Changing DEC");
	DEBUG_PRINTLN(dec_deg);
	dec_deg = (secs / 3600.0 + mins / 60.0 + deg) * multi;
	DEBUG_PRINTLN(dec_deg);
	isPositiveDeclination = multi > 0;

	// Store the new target declination.
	futurePosition.declination = dec_deg;
	
	// If there is a target select pin we need to reset the selected position to "none"
	#ifdef TARGET_SELECT_PIN
		selectedDebugTargetIndex = -1;
	#endif
}


/**
 * This gets called whenever a complete command was received.
 * It parses the received characters and calls the required functions
 * TODO Break this up into small functions so that the code stays maintainable
 */
bool parseCommands(Mount &telescope, Observer& observer) {
	if (newData == true) {
		if (receivedChars[0] == 'G' && receivedChars[1] == 'R') {
			// GR: Get Right Ascension
			getRightAscension(telescope);
		} else if (receivedChars[0] == 'G' && receivedChars[1] == 'D') {
			// GD: Get Declination
			getDeclination(telescope);
		} else if (receivedChars[0] == 'Q') {
			// Quit the current move
			moveQuit(telescope);
		} else if (receivedChars[0] == 'M' && receivedChars[1] == 'S') {
			// MS: Move Start
			// The function returning true means that isAligned was set to true.
			if (moveStart(telescope)) {
				// Ignores the next movement update and treats it as a home command instead
				//telescope.ignoreUpdates();

				newData = false;
				// Aligning was just performed
				return true;
			}
			else {
				// The telescope should not ignore movement updates once homed
				telescope.ignoreUpdates(false);
			}
		} else if (receivedChars[0] == 'S' && receivedChars[1] == 'r') {
			// Set Right Ascension (in hours, minutes and seconds)
			//telescope.ignoreUpdates();
			setRightAscension(telescope);
			display_statusUpdate(telescope);
		} else if (receivedChars[0] == 'S' && receivedChars[1] == 'd') {
			// Set target Declination (in degrees, minutes and seconds)
			//telescope.ignoreUpdates();
			setDeclination(telescope);
			display_statusUpdate(telescope);
		}
		else if (receivedChars[0] == 'T' && receivedChars[1] == 'R' && receivedChars[2] == 'K') {
			// Enable / Disable tracking (= home off / on)
			if (receivedChars[3] == '1') {
				telescope.setHomed(true);
				Serial.println("Enabled tracking");

				newData = false;
				return true;
			}
			else if (receivedChars[3] == '0') {
				telescope.setHomed(false);
				Serial.println("Disabled tracking");
			}
		}
		else if (receivedChars[0] == 'S' && receivedChars[1] == 'T' && receivedChars[2] == 'P') {
			// Enable / Disable stepper motors
			if (receivedChars[3] == '1') {
				digitalWrite(ALT_ENABLE_PIN, LOW);
				digitalWrite(AZ_ENABLE_PIN, LOW);
				Serial.println("Enabled stepper motors. Send :STP0# to disable them");
			}
			else if (receivedChars[3] == '0') {
				digitalWrite(ALT_ENABLE_PIN, HIGH);
				digitalWrite(AZ_ENABLE_PIN, HIGH);
				Serial.println("Disabled stepper motors. Send :STP1# to re-enable them");
			}
		}
		else if (receivedChars[0] == 'D' && receivedChars[1] == 'B'
				&& receivedChars[2] == 'G') {
			// DEBUG messages
			if (receivedChars[3] == 'M') {
				if (receivedChars[4] == 'I' || receivedChars[4] == 'D') {
					// Increase or Decrease Azimuth or Declination
					int add = receivedChars[4] == 'I' ? 1 : -1;
					if (receivedChars[5] == 'A') {
						// DBGMAXXX Move Azimuth to XXX
						ra_deg += add;
						RaDecPosition pos = telescope.getTarget();
						pos.rightAscension += add;
						telescope.setTarget(pos);
						Serial.println(
								add > 0 ?
										"Add 1 deg ascension" :
										"Sub 1 deg ascension");
						display_statusUpdate(telescope);
					} else if (receivedChars[5] == 'D') {
						// DBGMD[+/-]XX Move Declination to +/-XX
						dec_deg += add;
						RaDecPosition pos = telescope.getTarget();
						pos.declination += add;
						telescope.setTarget(pos);
						Serial.println(
								add > 0 ?
										"Add 1 deg declination" :
										"Sub 1 deg declination");
						display_statusUpdate(telescope);
					}
				} else {
					// Debug move to position stored in debugPositions[targetIndex]
					int targetIndex = char_to_int(receivedChars[4]);
					if (targetIndex > maxDebugPos) {
						Serial.println("Invalid index");
					} else {
						Serial.println();
						Serial.println("-----------------------------------------");
						Serial.println("Moving telescope to new target");

						Serial.print("Name\t");
						Serial.println(positionNames[targetIndex]);
						Serial.print("Ra\t");
						Serial.print(debugPositions[targetIndex][0]);
						Serial.println("�");
						Serial.print("Dec\t");
						Serial.print(debugPositions[targetIndex][1]);
						Serial.println("�");


						Serial.println("-----------------------------------------");

						ra_deg = debugPositions[targetIndex][0];
						dec_deg = debugPositions[targetIndex][1];
						RaDecPosition newPos = { debugPositions[targetIndex][0], debugPositions[targetIndex][1] };
						Serial.print("Scope msg: ");
						telescope.setTarget(newPos);
						Serial.println();
						Serial.println();
						display_statusUpdate(telescope);

						// If there is a target select button we need to store the selected position in 
						#ifdef TARGET_SELECT_PIN
							selectedDebugTargetIndex = targetIndex;
						#endif
					}
				}
			} else if (receivedChars[3] == 'H') {
				// Home the scope
				telescope.setHomed(true);
			} else if (receivedChars[3] == 'D' && receivedChars[4] == 'M') {
				// Disable Motors and Pause for X seconds
				digitalWrite(ALT_ENABLE_PIN, HIGH);
				digitalWrite(AZ_ENABLE_PIN, HIGH);
				Serial.print("Disabling motors for: ");
				const int disable_seconds = multi_char_to_int(receivedChars[5], receivedChars[6]);
				Serial.println(disable_seconds);
				delay(disable_seconds * 1000);
				Serial.println("Continuing");
				digitalWrite(ALT_ENABLE_PIN, LOW);
				digitalWrite(AZ_ENABLE_PIN, LOW);
			} else if (receivedChars[3] == 'G' && receivedChars[4] == 'P' && receivedChars[5] == 'S') {
				// Observer/Gps Debug info
				observer.printDebugInfo();
			}
			#ifdef SERIAL_DISPLAY_ENABLED
				else if (receivedChars[3] == 'D' && receivedChars[4] == 'S' && receivedChars[5] == 'P') {
					// Send the "Status: Online" command to the display
					Serial.println("Sending status update (online) to display");
					display_statusUpdate(telescope);
					Serial.println("Done...");
				}
			#endif
		} else if (receivedChars[0] == 'H' && receivedChars[1] == 'L' && receivedChars[2] == 'P') {
			printHelp();
		} else {
			Serial.println("ERROR: Unknown command");
			Serial.println(receivedChars);
		}
		newData = false;
	}
	return false;
}

bool handleSerialCommunication(Mount &telescope, Observer &observer) {
	// Receives the next command character, if available
	receiveCommandChar();

	// Once a complete command has been received, this parses and runs the command
	// Returns true, if the received command triggered a change from Mode::ALIGNING to Mode::TRACKING
	bool switchedToTracking = parseCommands(telescope, observer);

	// If we have a target select button it gets handled here
	#ifdef TARGET_SELECT_PIN
		if (!targetSelectButtonPressed && digitalRead(TARGET_SELECT_PIN) == LOW) {
			// Only one button click per 0.5 seconds is valid
			if (timeOfLastTargetChange + 500 < millis()) {
				timeOfLastTargetChange = millis();

				// Track the button state
				targetSelectButtonPressed = true;

				// Increase the selected target by 1 and wrap around to 0 if the result is larger than the max index
				selectedDebugTargetIndex = (selectedDebugTargetIndex + 1) % (maxDebugPos + 1);

				ra_deg = debugPositions[selectedDebugTargetIndex][0];
				dec_deg = debugPositions[selectedDebugTargetIndex][1];
				RaDecPosition newPos = { debugPositions[selectedDebugTargetIndex][0], debugPositions[selectedDebugTargetIndex][1] };
				telescope.setTarget(newPos);

				// Print confirmation and buzz
				DEBUG_PRINTLN("Switching target to " + String(selectedDebugTargetIndex));
				#ifdef BUZZER_PIN
					digitalWrite(BUZZER_PIN, HIGH);
					delay(100);
					digitalWrite(BUZZER_PIN, LOW);
				#endif
			}
		}

		if (targetSelectButtonPressed && digitalRead(TARGET_SELECT_PIN) == HIGH) {
			// Button released
			targetSelectButtonPressed = false;
		}
	#endif

	return switchedToTracking;
}