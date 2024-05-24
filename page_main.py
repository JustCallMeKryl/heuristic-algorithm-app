import streamlit as st
from streamlit_option_menu import option_menu

# Импорт страниц
import page_Plotnikov_Zverev_criterions
import page_Goldberg
import page_Plotnikov_Zverev_criterions_plus_Goldberg
import page_analysis


# Главная функция для навигации между страницами
def main():
    st.sidebar.title("Навигация")
    # Создание бокового меню
    with st.sidebar:
        selected = option_menu(
            menu_title="",
            options=["Главная", "Алгоритм Плотникова-Зверева", "Алгоритм Голдберга",
                     "Алгоритм Плотникова-Зверева + Голдберг", "Анализ"],
            icons=["house", "graph-up", "book", "bar-chart", "search"],
            menu_icon="cast",
            default_index=0,
        )

    # Логика навигации
    if selected == "Главная":
        st.title('Главная страница')
    elif selected == "Алгоритм Плотникова-Зверева":
        page_Plotnikov_Zverev_criterions.app()
    elif selected == "Алгоритм Голдберга":
        page_Goldberg.app()
    elif selected == "Алгоритм Плотникова-Зверева + Голдберг":
        page_Plotnikov_Zverev_criterions_plus_Goldberg.app()
    elif selected == "Анализ":
        page_analysis.app()


# Вызов главной функции
if __name__ == "__main__":
    main()
