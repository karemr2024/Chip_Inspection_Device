const int LED1=8;
const int LED2=9;
const int LED3=10;

int VAL=0;

void setup() 
   { 
      Serial.begin(9600); 
      pinMode(LED1, OUTPUT);
      pinMode(LED2, OUTPUT);
      pinMode(LED3, OUTPUT);
      digitalWrite (LED1, LOW);
      digitalWrite (LED2, LOW);
      digitalWrite (LED3, LOW);
   }
 
void loop() 
   {
     while (Serial.available())
        {
           VAL = Serial.read();
        }    
     if (VAL == '0'){
        digitalWrite (LED1, LOW);
        digitalWrite (LED2, LOW);
        digitalWrite (LED3, LOW);
     }  
     else if (VAL == '1'){
        digitalWrite (LED1, HIGH);
        digitalWrite (LED2, LOW);
        digitalWrite (LED3, LOW);
     }
     else if (VAL == '2'){
        digitalWrite (LED1, LOW);
        digitalWrite (LED2, HIGH);
        digitalWrite (LED3, LOW);
     } 
     else if (VAL == '3'){
        digitalWrite (LED1, LOW);
        digitalWrite (LED2, LOW);
        digitalWrite (LED3, HIGH);
     } 
   }
