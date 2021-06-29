/* Create a WiFi access point and provide a web server on it. */

#include <ESP8266WiFi.h>
#include <WiFiClient.h>
#include <ESP8266WebServer.h>
#include <DHT.h>
#include <DHT_U.h>
#include <Adafruit_Sensor.h>
#include <WebSocketsServer.h>

#define DHT_PIN 0
#define DHTTYPE DHT11

DHT dht(DHT_PIN, DHTTYPE);

float t =0.0;
float h = 0.0;

//SimpleDHT11 dht11(2); USAR GPIO0

#define TRIGGER_PIN  1  // Arduino pin tied to trigger pin on the ultrasonic sensor.
#define ECHO_PIN     2  // Arduino pin tied to echo pin on the ultrasonic sensor.
int cnt=0;

#ifndef APSSID
#define APSSID "ESPap"
#define APPSK  "thereisnospoon"
#endif

/* Set these to your desired credentials. */
const char *ssid = APSSID;
const char *password = APPSK;

ESP8266WebServer server(80);

/* Just a little test message.  Go to http://192.168.4.1 in a web browser
   connected to this access point to see it.
*/

void handleRoot() {
      String out;

      t = dht.readTemperature();
      h = dht.readHumidity();
      digitalWrite(TRIGGER_PIN, HIGH);  // SONAR trigger
      delayMicroseconds( 10 );       // wait 10 us
      digitalWrite(TRIGGER_PIN, LOW);  
     
      /*  echo time, t_eco (us) */
      unsigned long t_eco = pulseIn(ECHO_PIN, HIGH);
      cnt = cnt +1;
      out = String(cnt) + ", " + String(t_eco) + "," + String(t) + "," + String(h);
      
      server.send(200, "text/html", out);
}

void setup() {
  delay(1000);

  //GPIO1 (TX) swap the pin to a GPIO (TX -> FUNCTION_0)
  pinMode(1, FUNCTION_3);
  //GPIO3 (RX) swap the pin to a GPIO (RX -> FUNCTION_0)
  pinMode(3, FUNCTION_3);
  
  pinMode(TRIGGER_PIN, OUTPUT);
  pinMode(ECHO_PIN, INPUT);

  //dht.begin();
  
  /* You can remove the password parameter if you want the AP to be open. */
  WiFi.softAP(ssid, password);

  IPAddress myIP = WiFi.softAPIP();
  server.on("/", handleRoot);
  server.begin();
}

void loop() {
  server.handleClient();
}
