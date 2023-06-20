// Submitted by
// MD ABU OBAIDA ZISHAN
// ID:18201214
// Please note that this code was developed in an ubuntu machine
// might not work in windows

#include <iostream>
#include <fstream>
#include <unistd.h>
using namespace std;
unsigned int ONE_SECOND = 1000000;
bool led_state = false;

// please configure this port from arduino ide (IDE>TOOLS>PORT)
void sendData(string command)
{
  std::string port = "/dev/ttyUSB0";
  std::ofstream ard(port); // open the device
  usleep(ONE_SECOND);
  if (ard)
  {
    ard << command << '\n'; // send the command
  }
  else
    cout << "error"
         << "\n";
}
int main()
{
  while (1)
  {
    sendData("H");
    usleep(ONE_SECOND);
  }
  return 0;
}