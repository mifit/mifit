import Qt 4.7

Item {
    id: messageArea
    width: text.width
    height: text.height
    x: (root.width - width)/2
    anchors.leftMargin: 5
    anchors.bottom: root.bottom
    anchors.bottomMargin: 10

    Text {
        id: text
        color: "white"
        text: message
    }
}
