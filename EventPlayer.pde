
String Port1 = "COM19";      //Change these to your port names
String Port2 = "COM10";      //They can be found in the console
String Port3 = "COM18";      //when you run the file
String ArduinoPort = "COM24";
boolean EnablePorts = true;

int     LongestEventLength = 6000;    //Just make this big enough to stop getting index out of bounds errors
float   NanoSecondsPerFrame  = 50;
int     LedOnTime = 60;
float   FramesBetweenMuons = 10;
int     WaitAfterEvent = 75;
int     WaitBeforeEvent = 75;

int     CurrentFrame = 0;
int     testPix = 0;
int     MaxEvents = 200;
float   MinFreq = 100;
float   MaxFreq = 6000;
int     WaitTimer = -1;
float   MaxNoiseVolume = .05;

import processing.serial.*;
import java.awt.Rectangle;
import oscP5.*;
import netP5.*;

OscP5 oscP5;
NetAddress superColliderLocation;
NetAddress tabletLocation = new NetAddress("192.168.43.74", 54321);
boolean tabconnected;
boolean scconnected;
int time;
PImage img;
int dataw = 60, datah = 128;
float gamma = 1.7;

int numPorts=0;  // the number of serial ports in use
int maxPorts=24; // maximum number of serial ports

Serial[] ledSerial = new Serial[maxPorts];     // each port's actual Serial port
Rectangle[] ledArea = new Rectangle[maxPorts]; // the area of the movie each port gets, in % (0-100)
boolean[] ledLayout = new boolean[maxPorts];   // layout of rows, true = even is left->right
PImage[] ledImage = new PImage[maxPorts];      // image sent to each port
int[] gammatable = new int[256];
int errorCount=0;
float framerate=0;
boolean StopBackground = false;

int     CurrentPixels[][] = new int[datah * dataw][3];
float   EventData[][][] = new float[MaxEvents][LongestEventLength][9];
boolean EventsPlaying[] = new boolean[MaxEvents];
int     EventLengths[] = new int[MaxEvents];
int     DeathTimers[] = new int[dataw * datah];
int     DomCount = 5160;
float   Doms[][] = new float[DomCount][4];
float   EndFrames[] = new float[MaxEvents];
float   IceProperties[][] = new float[DomCount][3];
float   IdealMuons[][] = new float[14][8];
float   IdealElectrons[][] = new float[14][5];
int     TriggerTimer = -1;
int     EventToPlay = 0;
int     EventNums[][] = {
  { 2, 4, 7,12,17,22,27,-1,-1,-1,-1,-1,-1,-1},  //real  muon
  {-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2},  //ideal e-
  { 0, 1, 3, 5, 6, 8, 9,10,11,13,14,15,16,18},  //real  e-
  {-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3},  //ideal tau
  {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},  //real  tau
  {-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4}   //ideal muon
};

Serial ArduinoSerial;

void setup() {
  
  String RawFile[][] = new String[1][];
  RawFile[0] = loadStrings("LedGeometry.txt");
  for (int i = 0; i < DomCount; i++) {
    Doms[i] = float(split(RawFile[0][i], ' '));
  }
  /*
  RawFile[0] = loadStrings("LedIceProperties.txt");
  for (int i = 0; i < DomCount; i++) {
    IceProperties[i] = float(split(RawFile[0][i], ' '));
  }
  */
  for (int i = 0; i < MaxEvents; i++) {
    EventsPlaying[i] = false;
  }
  /*
  for (int i = 0; i < 14; i++) {
    float a = random(0, PI * 2);
    float b = random(-1, 1);
    IdealMuons[i][0] = random(-500,500);
    IdealMuons[i][1] = random(-500,500);
    IdealMuons[i][2] = random(-500,500);
    IdealMuons[i][3] = IdealMuons[i][0] + (2000 * cos(a) * sqrt(1 - pow(b,2)));
    IdealMuons[i][4] = IdealMuons[i][1] + (2000 * sin(a) * sqrt(1 - pow(b,2)));
    IdealMuons[i][5] = IdealMuons[i][2] + (2000 * b);
    IdealMuons[i][6] = random(30) + 120;
    IdealMuons[i][7] = 1.0;
    IdealElectrons[i][0] = random(-500,500);
    IdealElectrons[i][1] = random(-500,500);
    IdealElectrons[i][2] = random(-500,500);
    IdealElectrons[i][3] = random(100) + 250;
    IdealElectrons[i][4] = 1.0;
  }
  */
  for (int i = 0; i < datah * dataw; i++) {    //start with all strings off
    CurrentPixels[i][0] = 0;
    CurrentPixels[i][1] = 0;
    CurrentPixels[i][2] = 0;
    DeathTimers[i] = -1;
  }
  
  
  String[] list = Serial.list();
  
  delay(20);
  
  println("Serial Ports List:");
  println(list);
  if (EnablePorts) {
    ArduinoSerial = new Serial(this, ArduinoPort, 9600);
    serialConfigure(Port1);
    serialConfigure(Port2);
    serialConfigure(Port3);
    if (errorCount > 0) exit();
  }
  
  for (int i=0; i < 256; i++) {
    gammatable[i] = (int)(pow((float)i / 255.0, gamma) * 255.0 + 0.5);
  }
  
  size(480, 400);
  
  oscP5 = new OscP5(this, 12345);
  scconnected = false;
  time = 0;
  
  tabconnected = false;
  OscMessage myMessage = new OscMessage("/hello");
  myMessage.add("192.168.43.125"); // the laptop IP
  oscP5.send(myMessage, tabletLocation);
}

void update() {
  if (EnablePorts) {
    if (ArduinoSerial.available() > 0) {
      int Type = ArduinoSerial.read();
      for (int i = 0; i < MaxEvents; i++) {
        EventsPlaying[i] = false;
      }
      if (Type == '1') {
        TriggerTimer = WaitBeforeEvent;
        StopBackground = true;
        EventToPlay = int(floor(random(9)));
      }
      else if (Type == '2') {
        TriggerTimer = WaitBeforeEvent;
        StopBackground = true;
        EventToPlay = int(floor(random(9))) + 9;
      }
      else if (Type == '3') {
        TriggerTimer = WaitBeforeEvent;
        StopBackground = true;
        EventToPlay = int(floor(random(9))) + 18;
      }
      else println("Trigger unknown type: " + Type);
      while (ArduinoSerial.available() > 0) {
        ArduinoSerial.read();
      }
    }
  }
  for (int i = 0; i < MaxEvents; i++) {
    if (EventsPlaying[i]) {
      for (int j = 0; j < EventLengths[i]; j++) {
        if (EventData[i][j][0] == CurrentFrame) {
          CurrentPixels[int(EventData[i][j][1])][0] = int(EventData[i][j][2]);
          CurrentPixels[int(EventData[i][j][1])][1] = int(EventData[i][j][3]);
          CurrentPixels[int(EventData[i][j][1])][2] = int(EventData[i][j][4]);
          DeathTimers[int(EventData[i][j][1])] = LedOnTime;
          if (EventData[i][j][5] != 0) {
            TriggerGrain(EventData[i][j][5], EventData[i][j][6], EventData[i][j][7], EventData[i][j][8]);
          }
        }
      }
      if (CurrentFrame > EndFrames[i]) {
        EventsPlaying[i] = false;
      }
    }
  }
  CurrentFrame += 1;
  for (int i = 0; i < dataw * datah; i++) {
    if (DeathTimers[i] > 0) DeathTimers[i]--;
    else if (DeathTimers[i] == 0) {
      CurrentPixels[i][0] = 0;
      CurrentPixels[i][1] = 0;
      CurrentPixels[i][2] = 0;
      DeathTimers[i] = -1;
    }
  }
  if (StopBackground) {
    if (TriggerTimer > 0) {
      TriggerTimer--;
    }
    else if (TriggerTimer == 0) {
      TriggerTimer = -1;
      for (int i = 0; i < datah * dataw; i++) {
        CurrentPixels[i][0] = 0;
        CurrentPixels[i][1] = 0;
        CurrentPixels[i][2] = 0;
        DeathTimers[i] = -1;
      }
      int num = int(floor(random(28)));
      LoadEvent(num);
    }
    else if (!EventsPlaying[0]) {
      StopBackground = false;
      WaitTimer = WaitAfterEvent;
    }
  }
  if (WaitTimer > 0) WaitTimer--;
  else if (WaitTimer == 0) {
    WaitTimer = -1;
  }
  else if (!StopBackground && WaitTimer == -1) {
    int numTries = int(random(100));
    for (int i = 0; i < numTries; i++) {
      if (random(100) < 5) {
        int num = floor(random(DomCount));
        CurrentPixels[int(Doms[num][0])][0] = int(random(255));
        CurrentPixels[int(Doms[num][0])][1] = int(random(255));
        CurrentPixels[int(Doms[num][0])][2] = int(random(255));
        DeathTimers[int(Doms[num][0])] = LedOnTime;
        TriggerGrain((2 - pow(2, Doms[num][2] / 1100.0 + .5)) * (MaxFreq - MinFreq) + MinFreq,
                     Doms[num][1] / 550,
                     Doms[num][3] / 550,
                     float(CurrentPixels[int(Doms[num][0])][0] + CurrentPixels[int(Doms[num][0])][1] + CurrentPixels[int(Doms[num][0])][2]) / 765.0 * MaxNoiseVolume
                     );
      }
    }
  }
  
  img = createImage(dataw, datah, RGB);
  
  for (int i = 0; i < img.pixels.length; i++) {
    img.pixels[i] = color(CurrentPixels[i][0], CurrentPixels[i][1], CurrentPixels[i][2]);
  }
  
  framerate = 30.0;
  
  for (int i=0; i < numPorts; i++) {    
    // copy a portion of the movie's image to the LED image
    int xoffset = percentage(img.width, ledArea[i].x);
    int yoffset = percentage(img.height, ledArea[i].y);
    int xwidth =  percentage(img.width, ledArea[i].width);
    int yheight = percentage(img.height, ledArea[i].height);
    ledImage[i].copy(img, xoffset, yoffset, xwidth, yheight,
                     0, 0, ledImage[i].width, ledImage[i].height);
    // convert the LED image to raw data
    byte[] ledData =  new byte[(ledImage[i].width * ledImage[i].height * 3) + 3];
    image2data(ledImage[i], ledData, ledLayout[i]);
    if (i == 0) {
      ledData[0] = '*';  // first Teensy is the frame sync master
      int usec = (int)((1000000.0 / framerate) * 0.75);
      ledData[1] = (byte)(usec);   // request the frame sync pulse
      ledData[2] = (byte)(usec >> 8); // at 75% of the frame time
    } else {
      ledData[0] = '%';  // others sync to the master board
      ledData[1] = 0;
      ledData[2] = 0;
    }
    // send the raw data to the LEDs  :-)
    ledSerial[i].write(ledData);
  }
}

void LoadEvent(int EventNumber) {
  int EmptyEvent = -1;
  for (int i = 0; i < MaxEvents; i++) {
    if (!EventsPlaying[i]) {
      EmptyEvent = i;
      i = MaxEvents;
    }
  }
  if (EmptyEvent != -1) {
    String RawFile[][] = new String[1][];
    int FileLength;
    RawFile[0] = loadStrings("../LedData/Event" + str(EventNumber) + ".txt");
    EventLengths[EmptyEvent] = RawFile[0].length;
    for (int i = 0; i < EventLengths[EmptyEvent]; i++) {
      EventData[EmptyEvent][i] = float(split(RawFile[0][i], ' '));
      EventData[EmptyEvent][i][0] = CurrentFrame + floor(EventData[EmptyEvent][i][0] / NanoSecondsPerFrame);
    }
    EndFrames[EmptyEvent] = EventData[EmptyEvent][EventLengths[EmptyEvent] - 1][0];
    EventsPlaying[EmptyEvent] = true;
  }
  else {
    println("Tried to load event " + str(EventNumber) + ", but memory is full");
  }
}

void CreateMuon(float InitX, float InitY, float InitZ, float FinalX, float FinalY, float FinalZ, float Radius, float Volume) {
  int EmptyEvent = -1;
  for (int i = 0; i < MaxEvents; i++) {
    if (!EventsPlaying[i]) {
      EmptyEvent = i;
      i = MaxEvents;
    }
  }
  if (EmptyEvent != -1) {
    float SimEventData[][] = new float[LongestEventLength][9];
    float A[] = new float[3];
    float B[] = new float[3];
    float C[] = new float[3];
    float D[] = new float[3];
    float E[] = new float[3];
    float CurrentRadius;
    
    B[0] = FinalX - InitX;
    B[1] = FinalY - InitY;
    B[2] = FinalZ - InitZ;
    
    EventLengths[EmptyEvent] = 0;
    for (int i = 0; i < DomCount; i++) {
      A[0] = Doms[i][1] - InitX;
      A[1] = Doms[i][2] - InitY;
      A[2] = Doms[i][3] - InitZ;
      
      C[0] = Dot(A, B) / Dot(B, B) * B[0];
      C[1] = Dot(A, B) / Dot(B, B) * B[1];
      C[2] = Dot(A, B) / Dot(B, B) * B[2];
      
      D[0] = A[0] - C[0];
      D[1] = A[1] - C[1];
      D[2] = A[2] - C[2];
      
      CurrentRadius = Magnitude(D);
      
      if (CurrentRadius < Radius) {
        E[0] = Doms[i][1] - FinalX;
        E[1] = Doms[i][2] - FinalY;
        E[2] = Doms[i][3] - FinalZ;
        if (Dot(E, B) / (Magnitude(E) * Magnitude(B)) < 0) {
          if (Dot(A, B) / (Magnitude(A) * Magnitude(B)) > 0) {
            EventData[EmptyEvent][EventLengths[EmptyEvent]][0] = (Magnitude(C) / 0.3) + (Magnitude(D) / 0.225);
            EventData[EmptyEvent][EventLengths[EmptyEvent]][1] = Doms[i][0];
            EventData[EmptyEvent][EventLengths[EmptyEvent]][2] = 1 - CurrentRadius / Radius;
            EventData[EmptyEvent][EventLengths[EmptyEvent]][3] = 1 - CurrentRadius / Radius;
            EventData[EmptyEvent][EventLengths[EmptyEvent]][4] = 1 - CurrentRadius / Radius;
//          EventData[EmptyEvent][EventLengths[EmptyEvent]][5] = 440.0 * pow(2.0, (floor((1.0 - (Doms[i][2] / 1100.0 + 0.5)) * 88.0) - 48.0) / 12.0);
            EventData[EmptyEvent][EventLengths[EmptyEvent]][5] = (2 - pow(2, Doms[i][2] / 1100.0 + 0.5)) * (MaxFreq - MinFreq) + MinFreq;
            EventData[EmptyEvent][EventLengths[EmptyEvent]][6] = Doms[i][1] / 550.0;
            EventData[EmptyEvent][EventLengths[EmptyEvent]][7] = Doms[i][3] / 550.0;
            EventData[EmptyEvent][EventLengths[EmptyEvent]][8] = (1 - CurrentRadius / Radius) * Volume;
            EventLengths[EmptyEvent]++;
          }
        }
      }
    }
    if (EventLengths[EmptyEvent] > 0) {
      float MinT = EventData[EmptyEvent][0][0];
      float MaxT = EventData[EmptyEvent][0][0];
      for (int i = 1; i < EventLengths[EmptyEvent]; i++) {
        if (EventData[EmptyEvent][i][0] < MinT) {
          MinT = EventData[EmptyEvent][i][0];
        }
        if (EventData[EmptyEvent][i][0] > MaxT) {
          MaxT = EventData[EmptyEvent][i][0];
        }
      }
      for (int i = 0; i < EventLengths[EmptyEvent]; i++) {
        EventData[EmptyEvent][i][0] = EventData[EmptyEvent][i][0] - MinT;
      }
      MaxT = MaxT - MinT;
      EndFrames[EmptyEvent] = floor(MaxT / NanoSecondsPerFrame) + CurrentFrame;
      for (int i = 0; i < EventLengths[EmptyEvent]; i++) {
        if (EventData[EmptyEvent][i][0] / MaxT * 2 < 1) {
          EventData[EmptyEvent][i][2] = floor(EventData[EmptyEvent][i][2] * (255 - (255 * EventData[EmptyEvent][i][0] / MaxT * 2)));
          EventData[EmptyEvent][i][3] = floor(EventData[EmptyEvent][i][3] * (255 * EventData[EmptyEvent][i][0] / MaxT * 2));
          EventData[EmptyEvent][i][4] = 0;
        }
        else {
          EventData[EmptyEvent][i][2] = 0;
          EventData[EmptyEvent][i][3] = floor(EventData[EmptyEvent][i][3] * (255 - (255 * ((EventData[EmptyEvent][i][0] / MaxT * 2) - 1))));
          EventData[EmptyEvent][i][4] = floor(EventData[EmptyEvent][i][4] * (255 * ((EventData[EmptyEvent][i][0] / MaxT * 2) - 1)));
        }
      }
      for (int i = 0; i < EventLengths[EmptyEvent]; i++) {
        EventData[EmptyEvent][i][0] = floor(EventData[EmptyEvent][i][0] / NanoSecondsPerFrame) + CurrentFrame;
      }
      EventsPlaying[EmptyEvent] = true;
    }
  }
  else {
    println("Tried to generate muon, but memory is full");
  }
}

void CreateElectron(float X, float Y, float Z, float R, float Volume) {
  int EmptyEvent = -1;
  for (int i = 0; i < MaxEvents; i++) {
    if (!EventsPlaying[i]) {
      EmptyEvent = i;
      i = MaxEvents;
    }
  }
  if (EmptyEvent != -1) {
    EventLengths[EmptyEvent] = 0;
    float A[] = new float[3];
    float CurrentRadius;
    for (int i = 0; i < DomCount; i++) {
      A[0] = Doms[i][1] - X;
      A[1] = Doms[i][2] - Y;
      A[2] = Doms[i][3] - Z;
      CurrentRadius = Magnitude(A);
      if (CurrentRadius < R) {
        EventData[EmptyEvent][EventLengths[EmptyEvent]][0] = CurrentRadius / 0.225;
        EventData[EmptyEvent][EventLengths[EmptyEvent]][1] = Doms[i][0];
        EventData[EmptyEvent][EventLengths[EmptyEvent]][2] = 1 - CurrentRadius / R;
        EventData[EmptyEvent][EventLengths[EmptyEvent]][3] = 1 - CurrentRadius / R;
        EventData[EmptyEvent][EventLengths[EmptyEvent]][4] = 1 - CurrentRadius / R;
        EventData[EmptyEvent][EventLengths[EmptyEvent]][5] = (2 - pow(2, Doms[i][2] / 1100.0 + 0.5)) * (MaxFreq - MinFreq) + MinFreq;
        EventData[EmptyEvent][EventLengths[EmptyEvent]][6] = Doms[i][1] / 550.0;
        EventData[EmptyEvent][EventLengths[EmptyEvent]][7] = Doms[i][3] / 550.0;
        EventData[EmptyEvent][EventLengths[EmptyEvent]][8] = (1 - CurrentRadius / R) * Volume;
        EventLengths[EmptyEvent]++;
      }
    }
    if (EventLengths[EmptyEvent] > 0) {
      float MinT = EventData[EmptyEvent][0][0];
      float MaxT = EventData[EmptyEvent][0][0];
      for (int i = 1; i < EventLengths[EmptyEvent]; i++) {
        if (EventData[EmptyEvent][i][0] < MinT) {
          MinT = EventData[EmptyEvent][i][0];
        }
        if (EventData[EmptyEvent][i][0] > MaxT) {
          MaxT = EventData[EmptyEvent][i][0];
        }
      }
      for (int i = 0; i < EventLengths[EmptyEvent]; i++) {
        EventData[EmptyEvent][i][0] = EventData[EmptyEvent][i][0] - MinT;
      }
      MaxT = MaxT - MinT;
      EndFrames[EmptyEvent] = floor(MaxT / NanoSecondsPerFrame) + CurrentFrame;
      for (int i = 0; i < EventLengths[EmptyEvent]; i++) {
        if (EventData[EmptyEvent][i][0] / MaxT * 2 < 1) {
          EventData[EmptyEvent][i][2] = floor(EventData[EmptyEvent][i][2] * (255 - (255 * EventData[EmptyEvent][i][0] / MaxT * 2)));
          EventData[EmptyEvent][i][3] = floor(EventData[EmptyEvent][i][3] * (255 * EventData[EmptyEvent][i][0] / MaxT * 2));
          EventData[EmptyEvent][i][4] = 0;
        }
        else {
          EventData[EmptyEvent][i][2] = 0;
          EventData[EmptyEvent][i][3] = floor(EventData[EmptyEvent][i][3] * (255 - (255 * ((EventData[EmptyEvent][i][0] / MaxT * 2) - 1))));
          EventData[EmptyEvent][i][4] = floor(EventData[EmptyEvent][i][4] * (255 * ((EventData[EmptyEvent][i][0] / MaxT * 2) - 1)));
        }
      }
      for (int i = 0; i < EventLengths[EmptyEvent]; i++) {
        EventData[EmptyEvent][i][0] = floor(EventData[EmptyEvent][i][0] / NanoSecondsPerFrame) + CurrentFrame;
      }
      EventsPlaying[EmptyEvent] = true;
    }
  }
  else {
    println("Tried to generate electron, but memory is full");
  }
}

float Dot(float VectorA[], float VectorB[]) {
  return(VectorA[0] * VectorB[0] + VectorA[1] * VectorB[1] + VectorA[2] * VectorB[2]);
}

float Magnitude(float VectorA[]) {
  return(sqrt(pow(VectorA[0], 2) + pow(VectorA[1], 2) + pow(VectorA[2], 2)));
}

void quicksort(int EventNum, int leftstart, int rightstart) {
  float pivot = EventData[EventNum][leftstart][0];
  int leftedge = leftstart;
  int rightedge = rightstart - 1;
  while (leftedge < rightedge) {
    if (EventData[EventNum][leftedge][0] >= pivot) {
      if (EventData[EventNum][rightedge][0] < pivot) {
        swap(EventNum, leftedge, rightedge);
        leftedge++;
        rightedge--;
      }
      else rightedge--;
    }
    else leftedge++;
  }
  if (pivot <= EventData[EventNum][rightedge][0]) rightedge--;
  swap(EventNum, rightedge, rightstart);
  if (leftstart < rightedge - 1) {
    quicksort(EventNum, leftstart, rightedge - 1);
  }
  if (rightedge + 1 < rightstart) {
    quicksort(EventNum, rightedge + 1, rightstart);
  }
}

void swap(int EventNum, int val1, int val2) {
  for (int k = 0; k < 9; k++) {
    float temp = EventData[EventNum][val1][k];
    EventData[EventNum][val1][k] = EventData[EventNum][val2][k];
    EventData[EventNum][val2][k] = temp;
  }
}

void draw() {
  update();

  // show the original video
  image(img, 0,80);
 
  // then try to show what was most recently sent to the LEDs
  // by displaying all the images for each port.
  for (int i=0; i < numPorts; i++) {
    // compute the intended size of the entire LED array
    int xsize = percentageInverse(ledImage[i].width, ledArea[i].width);
    int ysize = percentageInverse(ledImage[i].height, ledArea[i].height);
    // computer this image's position within it
    int xloc =  percentage(xsize, ledArea[i].x);
    int yloc =  percentage(ysize, ledArea[i].y);
    // show what should appear on the LEDs
    image(ledImage[i], 240 - xsize / 2 + xloc, 10 + yloc);
  }
}

///////////////////////////////////////////////////////////
// image2data converts an image to OctoWS2811's raw data format.
// The number of vertical pixels in the image must be a multiple
// of 8.  The data array must be the proper size for the image.
void image2data(PImage image, byte[] data, boolean layout) {
  int offset = 3;
  int x, y, xbegin, xend, xinc, mask;
  int linesPerPin = image.height / 8;
  int pixel[] = new int[8];
  
  for (y = 0; y < linesPerPin; y++) {
    if ((y & 1) == (layout ? 0 : 1)) {
      // even numbered rows are left to right
      xbegin = 0;
      xend = image.width;
      xinc = 1;
    } else {
      // odd numbered rows are right to left
      xbegin = image.width - 1;
      xend = -1;
      xinc = -1;
    }
    for (x = xbegin; x != xend; x += xinc) {
      for (int i=0; i < 8; i++) {
        // fetch 8 pixels from the image, 1 for each pin
        pixel[i] = image.pixels[x + (y + linesPerPin * i) * image.width];
        pixel[i] = colorWiring(pixel[i]);
      }
      // convert 8 pixels to 24 bytes
      for (mask = 0x800000; mask != 0; mask >>= 1) {
        byte b = 0;
        for (int i=0; i < 8; i++) {
          if ((pixel[i] & mask) != 0) b |= (1 << i);
        }
        data[offset++] = b;
      }
    }
  } 
}

 ///////////////////////////////////////////////////////////
// translate the 24 bit color from RGB to the actual
// order used by the LED wiring.  GRB is the most common.
int colorWiring(int c) {
  int red = (c & 0xFF0000) >> 16;
  int green = (c & 0x00FF00) >> 8;
  int blue = (c & 0x0000FF);
  red = gammatable[red];
  green = gammatable[green];
  blue = gammatable[blue];
  return (green << 16) | (red << 8) | (blue); // GRB - most common wiring
}


// ask a Teensy board for its LED configuration, and set up the info for it.
void serialConfigure(String portName) {
  if (numPorts >= maxPorts) {
    println("too many serial ports, please increase maxPorts");
    errorCount++;
    return;
  }
  try {
    ledSerial[numPorts] = new Serial(this, portName);
    if (ledSerial[numPorts] == null) throw new NullPointerException();
    ledSerial[numPorts].write('?');
  } catch (Throwable e) {
    println("Serial port " + portName + " does not exist or is non-functional");
    errorCount++;
    return;
  }
  delay(50);
  String line = ledSerial[numPorts].readStringUntil(10);
  if (line == null) {
    println("Serial port " + portName + " is not responding.");
    println("Is it really a Teensy 3.0 running VideoDisplay?");
    errorCount++;
    return;
  }
  String param[] = line.split(",");
  if (param.length != 12) {
    println("Error: port " + portName + " did not respond to LED config query");
    errorCount++;
    return;
  }
  // only store the info and increase numPorts if Teensy responds properly
  ledImage[numPorts] = new PImage(Integer.parseInt(param[0]), Integer.parseInt(param[1]), RGB);
  ledArea[numPorts] = new Rectangle(Integer.parseInt(param[5]), Integer.parseInt(param[6]),
                     Integer.parseInt(param[7]), Integer.parseInt(param[8]));
  ledLayout[numPorts] = (Integer.parseInt(param[5]) == 0);
  numPorts++;
}


// scale a number by a percentage, from 0 to 100
int percentage(int num, int percent) {
  double mult = percentageFloat(percent);
  double output = num * mult;
  return (int)output;
}

// scale a number by the inverse of a percentage, from 0 to 100
int percentageInverse(int num, int percent) {
  double div = percentageFloat(percent);
  double output = num / div;
  return (int)output;
}

// convert an integer from 0 to 100 to a float percentage
// from 0.0 to 1.0.  Special cases for 1/3, 1/6, 1/7, etc
// are handled automatically to fix integer rounding.
double percentageFloat(int percent) {
  if (percent == 33) return 1.0 / 3.0;
  if (percent == 17) return 1.0 / 6.0;
  if (percent == 14) return 1.0 / 7.0;
  if (percent == 13) return 1.0 / 8.0;
  if (percent == 11) return 1.0 / 9.0;
  if (percent ==  9) return 1.0 / 11.0;
  if (percent ==  8) return 1.0 / 12.0;
  return (double)percent / 100.0;
}

void eventPlay(int eventType, int eventNum) {
  println("/event " + str(eventType) + " " + str(eventNum));
  OscMessage myMessage = new OscMessage("/eventstop");
  myMessage.add(eventType);
  myMessage.add(eventNum);
  oscP5.send(myMessage, tabletLocation);
  if (EventNums[eventType][eventNum] >= 0) {
    LoadEvent(EventNums[eventType][eventNum]);
  }
  else if (EventNums[eventType][eventNum] == -2) {
    CreateElectron(random(-500,500),random(-500,500),random(-500,500),random(50)+200,1.0);
  }
  else if (EventNums[eventType][eventNum] == -4) {
    CreateMuon(
      IdealMuons[eventNum][0],
      IdealMuons[eventNum][1],
      IdealMuons[eventNum][2],
      IdealMuons[eventNum][3],
      IdealMuons[eventNum][4],
      IdealMuons[eventNum][5],
      IdealMuons[eventNum][6],
      IdealMuons[eventNum][7]
    );
  }
  else if (EventNums[eventType][eventNum] == -3) {
    CreateElectron(
      IdealElectrons[eventNum][0],
      IdealElectrons[eventNum][1],
      IdealElectrons[eventNum][2],
      IdealElectrons[eventNum][3],
      IdealElectrons[eventNum][4]
    );
  }
}

/////OSC MESSAGE HANDLER/////////////////////////////////////////////////////////////
/// incoming osc message are forwarded to the oscEvent method.
void oscEvent(OscMessage theOscMessage) {
// print the address pattern and the typetag of the received OscMessage
//  print("### received an osc message.");
//  print(" addrpattern: "+theOscMessage.addrPattern());
//  println(" typetag: "+theOscMessage.typetag());
//  println(" typetag: "+theOscMessage.netAddress().port());
  
  ///// event
  if (theOscMessage.checkAddrPattern("/event")) {
    if (theOscMessage.checkTypetag("ii")) {
      // get the event number
      int eventtype = (int)theOscMessage.get(0).intValue(); //comes from control as int
      int eventnum = (int)theOscMessage.get(1).intValue(); //comes from control as int
      
      //trigger the event
      eventPlay(eventtype, eventnum);
    }
  ///// event count
  } else if (theOscMessage.checkAddrPattern("/speed")) {
    float msgVal = theOscMessage.get(0).floatValue();
    NanoSecondsPerFrame = msgVal;
  
  ///// grainsize
  } else if (theOscMessage.checkAddrPattern("/grainsize")) {
    if (scconnected) {
      float msgVal = theOscMessage.get(0).floatValue();
      OscMessage myMessage = new OscMessage("/grainsize");
      myMessage.add(msgVal);
      oscP5.send(myMessage, superColliderLocation);
      
      println("/grainsize " + msgVal);
    } else {
      println("/grainsize, sc not connected");
    }
  
  ///// volume
  } else if (theOscMessage.checkAddrPattern("/volume")) {

    if (scconnected) {
      float msgVal = theOscMessage.get(0).floatValue();
      OscMessage myMessage = new OscMessage("/volume");
      myMessage.add(msgVal);
      oscP5.send(myMessage, superColliderLocation);
      println("/volume " + msgVal);
    } else {
      println("/volume, sc not connected");
    }
  
  //scip
  } else if (theOscMessage.checkAddrPattern("/scip")) {
  
  
  /* superColliderLocation is a NetAddress. a NetAddress takes 2 parameters,
   * an ip address and a port number. superColliderLocation is used as parameter in
   * oscP5.send() when sending osc packets to another computer, device, application. 
   * For usage see above.
   */
    superColliderLocation = new NetAddress(theOscMessage.netAddress().address(), 57120);//the address of SuperCollider
    
    OscMessage myMessage = new OscMessage("/ack");
    oscP5.send(myMessage, superColliderLocation);
    
    println("connected");
    scconnected = true;
  } else if (theOscMessage.checkAddrPattern("/tabletip")) {
   if (theOscMessage.checkTypetag("i")) {
      println("/tabletip " + (int)theOscMessage.get(0).intValue());
      
      // populate the event players
      OscMessage myMessage1 = new OscMessage("/eventcount");
      myMessage1.add(0);
      myMessage1.add(7);//real muon
      oscP5.send(myMessage1, tabletLocation);
      
      OscMessage myMessage2 = new OscMessage("/eventcount");
      myMessage2.add(1);
      myMessage2.add(14);//ideal e-
      oscP5.send(myMessage2, tabletLocation);
      
      OscMessage myMessage3 = new OscMessage("/eventcount");
      myMessage3.add(2);
      myMessage3.add(14);//real e-
      oscP5.send(myMessage3, tabletLocation);
      
      OscMessage myMessage4 = new OscMessage("/eventcount");
      myMessage4.add(3);
      myMessage4.add(0);//ideal tau
      oscP5.send(myMessage4, tabletLocation);
      
      OscMessage myMessage5 = new OscMessage("/eventcount");
      myMessage5.add(4);
      myMessage5.add(0);//real tau
      oscP5.send(myMessage5, tabletLocation);
      
      OscMessage myMessage6 = new OscMessage("/eventcount");
      myMessage6.add(5);
      myMessage6.add(14);//ideal muon
      oscP5.send(myMessage6, tabletLocation);
    } else {
      println("error: event type tag mismatch - expecting type i, received "
      +theOscMessage.typetag());
    }
  }
}

void TriggerGrain (float freq, float x, float y, float amp) {
 if (scconnected) {
    OscMessage myMessage = new OscMessage("/grain");
    myMessage.add(freq);
    myMessage.add(x);
    myMessage.add(y);
    myMessage.add(amp);
    oscP5.send(myMessage, superColliderLocation);
  }
}

