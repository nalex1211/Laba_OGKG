module org.example.laba_ogkg {
    requires javafx.controls;
    requires javafx.fxml;
    requires commons.math3;


    opens org.example.laba_ogkg to javafx.fxml;
    exports org.example.laba_ogkg;
}