#include<iostream>
#include <windows.h>
#include<cstdio>

using namespace std;

int main()
{
    HWND h = ::FindWindow(NULL, "Yuri's Revenge");
    DWORD processid;
    GetWindowThreadProcessId(h, &processid);
    HANDLE hprocess = 0;

    hprocess = OpenProcess(PROCESS_ALL_ACCESS, FALSE, processid);

	if (hprocess == 0)
    {
		cout << "Opened False!" << endl;
		return 1;
	} 
	else
	{
		//修改钱的值
		DWORD money = 123;
		//00a82cb4+30c 基址+偏移
        DWORD addr = 0x00a82cb4;

        ReadProcessMemory(hprocess,(LPVOID)(addr), &addr, 4, NULL);
		//偏移量
        addr += 0x30C;
        WriteProcessMemory(hprocess, (LPVOID)addr, &money, 4, 0);

        cout << "Successed!" << endl;

		return 0;
	}
    return 0;
}
