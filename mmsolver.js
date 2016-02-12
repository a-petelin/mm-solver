/*
1. Прямо
	а. Переложить 
	б. Построение решения гиперболических уравнений по вычисленным потокам
	в. Рассчёт g и f-критериев 
	г. Проставить сеточный флаг "есть точки по критериям"
	д. Решение эллиптических уравнений (предыдущее решение запоминается)
	е. Рассчёт потоков для гиберболических уравнений 
2. Вниз
	а. Помещение точек согласно g,f-критериям в буфер, ели они есть
	б. Помещение граничных точек в буфер (граничная точка сетки -- это такая точка не удолетворяющая g,f-критериям, но имеющая таковую в соседях по Муру на расстоянии 1)
	*. При помещении точек вместо решений эллиптических уравнений забирается усреднение этого решения по граням.
	в. Разбить каждую точку в буфере на куб из k^d точек.
	г. 
	д. Увеличить номер текущей сетки на 1
	е. Вставить точки из буфера в сетку
	ё. Разбить сетку на блоки
	ж. Для точек каждого блока остаить только те значения эллиптических уравнений, которые лежать на граничных для блока гранях.
3. Вверх
	а. 
	б. Поместить все неграничные точки сетки в буфер
	в. Каждый куб из k^d точек буфера собрать в одну точку
	г. Удалить точки согласно g,f-критериям 	
	д. Уменьшить номер текущей сетки на 1.
	е. Вставить точки из буфера
	В точках граничных с вставленными потоки в эти точки установить согласно значениям на соответствующих гранях вставленных точек
	
Изначально считается что был переход прямо
Если предыдущий переход был прямо, то 
-- Если "есть точки по критериям", то делаем переход вниз
-- Если 
*/