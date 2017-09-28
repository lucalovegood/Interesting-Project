#include<windows.h>
#include<cstdio>
#include<iostream>

using namespace std;

HWND hq=FindWindow(NULL, "冒险岛2 - 尬萌不限号");

int spaceBar = 32;

int main()
{

//    HWND hq=FindWindow(NULL, "a.txt - 记事本");

    if (hq != NULL)
    {
        cout << "找到对应窗口" << endl;
//        ShowWindow(hp, SW_RESTORE);
       //设置为活动窗体，防止被其他窗口挡住
//       SetForegroundWindow(hp);
        while(true)
        {
            SendMessage(hq, WM_KEYDOWN, spaceBar, 0);
            keybd_event(spaceBar, MapVirtualKey(spaceBar, 0), 0 ,0);
            Sleep(300);
            keybd_event(spaceBar, MapVirtualKey(spaceBar, 0), KEYEVENTF_KEYUP ,0);
            printf("a\n");
            Sleep(3000);
        }
    }
    else
    {
        cout << "没找到对应窗口" << endl;
    }

    return 0;
}
