# Параллельная обработка изображений (Бакалаврская выпускная работа)

Здесь Вы можете посмотреть мою бакалаврскую выпускную квалификационную работу, написанную на языке C++ с применением
объектно-ориентированного программирования(ООП).
Тема работы "Разработка алгоритма и программы для параллельной обработки изображений".
В качестве алгоритма для параллельной обработки данных/изображений используется Обход графа в ширину(Поиск в ширину).
В качестве технологии параллельных вычислений для ускорения процесса обработки данных используется MPI.
В качестве xml-парсера используется RapidXML, так как это самый быстрый парсер на сегодняшний день.
В программе реализованы 3 класса объектов: Графы, Вершины и Соседи. У каждого класса свои поля и методы.

Чтобы проверить работу программы самостоятельно, нужно:
1) Установить Visual Studio;
2) MPI нужно установить как показано в видео https://www.youtube.com/watch?v=dhxW4ZoZQdI&t=23s&ab_channel=LandLord
3) Запустить ConsoleApplication1.sln;
4) Выбрать конфигурацию решения Release;
5) Нажать Локальный отладчик Windows;
6) Поставить в появившемся окне галочки Разрешить связь в частных и общественных сетях;
7) Завершить отладку;
8) Запустить Total Commander(так удобнее), указать путь до папки ...\x64\Release\
9) В самом низу Total Commander'a в белой пустой строке написать cmd
10) В появившеся окне написать или скопировать вставить "mpiexec -n 16 ConsoleApplication1" без кавычек.
11) Зайдите в папку Output files, проверьте появившиеся файлы.